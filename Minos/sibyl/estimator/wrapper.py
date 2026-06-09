import torch
import torch.nn.functional as F
from pathlib import Path
from typing import Dict, Tuple

from uGUIDE.config_utils import create_config_uGUIDE, save_config_uGUIDE, load_config_uGUIDE
from uGUIDE.training import run_training
from uGUIDE.estimation import estimate_microstructure
from uGUIDE.embedded_net import get_embedded_net
from uGUIDE.normalization import load_normalizer

from sibyl.data.synthetic import PRIOR_BOUNDS

def create_and_train_estimator(
    theta_train: torch.Tensor,
    x_train: torch.Tensor,
    model_name: str = 'ivim_synthetic',
    folderpath: str = './results',
    epochs: int = 100,
    seed: int = 42
) -> Dict:
    """
    Creates uGUIDE config, trains the NPE, and returns the config.
    """
    folderpath_p = Path(folderpath)
    folderpath_p.mkdir(parents=True, exist_ok=True)
    
    config = create_config_uGUIDE(
        microstructure_model_name=model_name,
        size_x=x_train.shape[1],
        prior=PRIOR_BOUNDS,
        folderpath=folderpath_p,
        max_epochs=epochs,
        random_seed=seed,
        use_MLP=True,
        nf_features=6, # 2 * size_theta = 6
    )
    
    save_config_uGUIDE(config, savefile='config.pkl', folderpath=folderpath_p)
    
    print("Training uGUIDE NPE...")
    run_training(theta_train, x_train, config=config, plot_loss=True, load_state=False)
    print("Training complete.")
    
    return config

def get_posterior(x: torch.Tensor, config: Dict) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Wrapper around uGUIDE's estimate_microstructure.
    Returns MAP and uncertainty (IQR).
    """
    # estimate_microstructure returns: map, mask, degeneracy_mask, uncertainty, ambiguity
    # Since we need posterior statistics, map and uncertainty (which is IQR) are returned.
    map_est, mask, degeneracy_mask, uncertainty, ambiguity = estimate_microstructure(
        x, config, verbose=False, plot=False
    )
    return map_est, uncertainty

def get_embedding(x: torch.Tensor, config: Dict) -> torch.Tensor:
    """
    Extracts the summary embedding from the trained uGUIDE model.
    Crucially, includes the exact same nonlinearities applied internally by uGUIDE.
    """
    device = config['device']
    x = x.to(device)
    
    x_normalizer = load_normalizer(config['folderpath'] / config['x_normalizer_file']).to(device)
    embedded_net = get_embedded_net(
        input_dim=config['size_x'],
        output_dim=config['nf_features'],
        layer_1_dim=config['hidden_layers'][0],
        layer_2_dim=config['hidden_layers'][1],
        pretrained_state=config['folderpath'] / config['embedder_state_dict_file'],
        use_MLP=config['use_MLP'],
        device=device
    ).to(device)
    
    embedded_net.eval()
    
    with torch.inference_mode():
        x_norm = x_normalizer(x)
        embedding = embedded_net(x_norm)
        
        # Exact operations performed internally by uGUIDE before normalizing flow
        embedding = torch.tanh(embedding)
        embedding = F.layer_norm(embedding, embedding.shape[-1:])
        
    return embedding

def test_embedding_matches_uguide(x: torch.Tensor, config: Dict):
    """
    Smoke test to verify that our get_embedding function identically matches 
    what uGUIDE's estimation module does internally.
    """
    import sys
    from unittest.mock import patch
    
    # We will patch torch.distributions.ConditionalTransformedDistribution.condition 
    # to intercept the embedding passed to the flow inside estimate_microstructure.
    intercepted_embedding = [None]
    
    original_condition = torch.distributions.ConditionalTransformedDistribution.condition
    
    def mock_condition(self, context):
        intercepted_embedding[0] = context.clone()
        return original_condition(self, context)
    
    with patch('torch.distributions.ConditionalTransformedDistribution.condition', mock_condition):
        _ = estimate_microstructure(x[:2], config, verbose=False, plot=False)
        
    internal_embedding = intercepted_embedding[0]
    wrapper_embedding = get_embedding(x[:2], config)
    
    assert internal_embedding is not None, "Smoke test failed: could not intercept internal embedding."
    assert torch.allclose(internal_embedding, wrapper_embedding, atol=1e-6), \
        f"Smoke test failed: wrapper embedding does not match internal uGUIDE embedding! Max diff: {torch.max(torch.abs(internal_embedding - wrapper_embedding))}"
    
    print("Smoke test passed: extracted embedding exactly matches flow condition context.")
