import torch
import numpy as np

class ResidualConformalDetector:
    """
    Family 2 Detector: Signal-space OOD detector using prediction residuals.
    Fits the ID distribution of signal residuals and scores test signals.
    """
    def __init__(self):
        self.calib_residuals = None

    def fit(self, id_residuals: torch.Tensor):
        """
        Fit the detector on ID residuals.
        
        Parameters
        ----------
        id_residuals : torch.Tensor
            Residuals of ID samples (shape: [N] or [N, B]).
            Usually we use the norm of the residual vector for each sample.
        """
        if id_residuals.ndim > 1:
            # Use L2 norm of the residual vector
            self.calib_residuals = torch.norm(id_residuals, p=2, dim=1).cpu().numpy()
        else:
            self.calib_residuals = id_residuals.cpu().numpy()
            
        # Sort for conformal p-value computation
        self.calib_residuals.sort()

    def score(self, test_residuals: torch.Tensor) -> torch.Tensor:
        """
        Compute the OOD score (residual magnitude) for test samples.
        Higher score means more OOD.
        
        Parameters
        ----------
        test_residuals : torch.Tensor
            Residuals of test samples (shape: [N] or [N, B]).
            
        Returns
        -------
        torch.Tensor
            OOD scores (shape: [N]).
        """
        assert self.calib_residuals is not None, "Detector must be fitted before scoring."
        
        if test_residuals.ndim > 1:
            test_res_norm = torch.norm(test_residuals, p=2, dim=1)
        else:
            test_res_norm = test_residuals
            
        # The score itself can just be the residual magnitude. 
        # (For conformal, we could return 1 - p_value, but magnitude is monotonically equivalent)
        return test_res_norm
