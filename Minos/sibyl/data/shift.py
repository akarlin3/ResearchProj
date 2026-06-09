import torch

def apply_noise(clean_signal: torch.Tensor, snr: float, rician: bool = False, seed: int = None) -> torch.Tensor:
    """
    Apply Gaussian or Rician noise to a clean signal.

    Parameters
    ----------
    clean_signal : torch.Tensor
        Clean IVIM signal.
    snr : float
        Signal-to-noise ratio (defined relative to S0=1.0).
    rician : bool, default=False
        If True, apply Rician noise. Otherwise, apply Gaussian noise.
    seed : int, optional
        Random seed.

    Returns
    -------
    torch.Tensor
        Noisy signal.
    """
    if seed is not None:
        torch.manual_seed(seed)
        
    sigma = 1.0 / snr
    
    noise_r = torch.randn_like(clean_signal) * sigma
    
    if not rician:
        return clean_signal + noise_r
    else:
        noise_i = torch.randn_like(clean_signal) * sigma
        return torch.sqrt((clean_signal + noise_r)**2 + noise_i**2)

def apply_corruption(noisy_signal: torch.Tensor, drop_prob: float, attenuation: float = 0.5, seed: int = None) -> torch.Tensor:
    """
    Apply multiplicative dropout corruption (simulating severe motion or artifact).
    Randomly selects b-values and attenuates them.

    Parameters
    ----------
    noisy_signal : torch.Tensor
        The signal to corrupt (shape: [N, B]).
    drop_prob : float
        Probability of corrupting a specific b-value for a given sample.
    attenuation : float, default=0.5
        The attenuation factor applied to the corrupted b-values.
        E.g. 0.5 means the signal drops by 50%.
    seed : int, optional
        Random seed.

    Returns
    -------
    torch.Tensor
        Corrupted signal.
    """
    if seed is not None:
        torch.manual_seed(seed)
        
    mask = torch.rand_like(noisy_signal) < drop_prob
    corrupted_signal = noisy_signal.clone()
    corrupted_signal[mask] = corrupted_signal[mask] * attenuation
    
    return corrupted_signal

def generate_shifted_datasets(theta: torch.Tensor, clean_signal: torch.Tensor, seed: int = 42):
    """
    Generates dictionary of OOD datasets according to Tier 1 shift axes.
    
    Axes:
    - Axis A (Noise Model): Gaussian -> Rician mismatch (at SNR 50).
    - Axis B (SNR): Rician SNR 50 -> 30 -> 15 -> 5.
    - Axis C (Corruption): Multiplicative dropout on Gaussian SNR 50.
      Levels vary the drop probability (10%, 20%, 30%, 40%).
      
    Returns
    -------
    dict
        Dictionary mapping shift name to tuple (theta, shifted_signal).
    """
    # Fix seed for reproducibility
    torch.manual_seed(seed)
    
    datasets = {}
    
    # ID is Gaussian SNR 50 (for reference, already generated in synthetic.py but keeping here for completeness)
    
    # Axis A: Noise Model Mismatch (and also serves as base for Axis B)
    rician_50 = apply_noise(clean_signal, snr=50.0, rician=True)
    datasets['axis_A_rician_snr50'] = (theta, rician_50)
    
    # Axis B: SNR Degradation (Rician)
    datasets['axis_B_rician_snr30'] = (theta, apply_noise(clean_signal, snr=30.0, rician=True))
    datasets['axis_B_rician_snr15'] = (theta, apply_noise(clean_signal, snr=15.0, rician=True))
    datasets['axis_B_rician_snr5'] = (theta, apply_noise(clean_signal, snr=5.0, rician=True))
    
    # Axis C: Corruption (applied on Gaussian SNR 50 to isolate the corruption effect)
    gaussian_50 = apply_noise(clean_signal, snr=50.0, rician=False)
    datasets['axis_C_corrupt_10'] = (theta, apply_corruption(gaussian_50, drop_prob=0.10, attenuation=0.5))
    datasets['axis_C_corrupt_20'] = (theta, apply_corruption(gaussian_50, drop_prob=0.20, attenuation=0.5))
    datasets['axis_C_corrupt_30'] = (theta, apply_corruption(gaussian_50, drop_prob=0.30, attenuation=0.5))
    datasets['axis_C_corrupt_40'] = (theta, apply_corruption(gaussian_50, drop_prob=0.40, attenuation=0.5))
    
    return datasets
