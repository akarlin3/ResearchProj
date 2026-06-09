import torch
from typing import Dict, Tuple, List
from sibyl.forward_model.ivim import ivim_biexponential, DEFAULT_B_SCHEME

PRIOR_BOUNDS = {
    'f': [0.01, 0.20],
    'D': [0.5e-3, 2.5e-3],
    'Dstar': [10.0e-3, 50.0e-3]
}

# Ensure the keys are ordered consistently
PRIOR_KEYS = ['f', 'D', 'Dstar']

def sample_prior(n_samples: int, seed: int = None) -> torch.Tensor:
    """
    Sample parameters uniformly from the predefined breast IVIM prior.

    Parameters
    ----------
    n_samples : int
        Number of samples to generate.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    torch.Tensor
        Sampled parameter tensor (shape: [n_samples, 3]). 
        Order is [f, D, Dstar].
    """
    if seed is not None:
        torch.manual_seed(seed)
        
    samples = []
    for key in PRIOR_KEYS:
        low, high = PRIOR_BOUNDS[key]
        # Uniform sampling
        param_samples = torch.empty(n_samples).uniform_(low, high)
        samples.append(param_samples)
        
    return torch.stack(samples, dim=1)

def generate_id_dataset(n_samples: int, snr: float = 50.0, seed: int = None) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Generate In-Distribution (ID) synthetic dataset.
    This includes sampling from the prior, simulating the IVIM signal, 
    and adding Gaussian noise at the specified SNR.

    Parameters
    ----------
    n_samples : int
        Number of samples to generate.
    snr : float, default=50.0
        Signal-to-noise ratio for Gaussian noise.
    seed : int, optional
        Random seed.

    Returns
    -------
    theta : torch.Tensor
        Parameter tensor (shape: [n_samples, 3])
    x : torch.Tensor
        Noisy signal tensor (shape: [n_samples, B])
    """
    if seed is not None:
        torch.manual_seed(seed)
        
    theta = sample_prior(n_samples, seed=seed)
    f, D, Dstar = theta[:, 0], theta[:, 1], theta[:, 2]
    
    clean_signal = ivim_biexponential(DEFAULT_B_SCHEME, f, D, Dstar)
    
    # Add Gaussian noise
    sigma = 1.0 / snr
    noise = torch.randn_like(clean_signal) * sigma
    noisy_signal = clean_signal + noise
    
    return theta, noisy_signal
