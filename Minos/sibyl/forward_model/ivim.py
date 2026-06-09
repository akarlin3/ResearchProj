import torch
from typing import Optional

# Default dense synthetic 10-value b-scheme for Tier 1
DEFAULT_B_SCHEME = torch.tensor([0, 10, 20, 50, 100, 200, 400, 600, 800, 1000], dtype=torch.float32)

def ivim_biexponential(
    bvals: torch.Tensor,
    f: torch.Tensor,
    D: torch.Tensor,
    Dstar: torch.Tensor,
    S0: Optional[torch.Tensor] = None
) -> torch.Tensor:
    """
    IVIM biexponential forward model.

    Parameters
    ----------
    bvals : torch.Tensor
        1D tensor of b-values (shape: [B])
    f : torch.Tensor
        Perfusion fraction (shape: [N, 1] or [N])
    D : torch.Tensor
        Tissue diffusion coefficient (shape: [N, 1] or [N])
    Dstar : torch.Tensor
        Pseudo-diffusion coefficient (shape: [N, 1] or [N])
    S0 : torch.Tensor, optional
        Signal at b=0. If None, S0=1.0 is used. (shape: [N, 1] or [N])

    Returns
    -------
    torch.Tensor
        Simulated signal (shape: [N, B])
    """
    if f.ndim == 1:
        f = f.unsqueeze(-1)
    if D.ndim == 1:
        D = D.unsqueeze(-1)
    if Dstar.ndim == 1:
        Dstar = Dstar.unsqueeze(-1)
    if S0 is not None and S0.ndim == 1:
        S0 = S0.unsqueeze(-1)
        
    bvals = bvals.view(1, -1)  # shape [1, B]
    
    signal = (f * torch.exp(-bvals * (D + Dstar))) + ((1 - f) * torch.exp(-bvals * D))
    
    if S0 is not None:
        signal = S0 * signal
        
    return signal
