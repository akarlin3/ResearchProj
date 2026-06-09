import torch
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
from scipy.stats import spearmanr
from typing import Dict, Tuple

def compute_detection_metrics(id_scores: np.ndarray, ood_scores: np.ndarray) -> Dict[str, float]:
    """
    Compute AUROC, AUPRC, and FPR@95%TPR for OOD detection.
    
    Parameters
    ----------
    id_scores : np.ndarray
        Scores for ID samples (lower is more ID).
    ood_scores : np.ndarray
        Scores for OOD samples (higher is more OOD).
        
    Returns
    -------
    Dict[str, float]
        Dictionary containing AUROC, AUPRC, and FPR95.
    """
    y_true = np.concatenate([np.zeros(len(id_scores)), np.ones(len(ood_scores))])
    y_scores = np.concatenate([id_scores, ood_scores])
    
    auroc = roc_auc_score(y_true, y_scores)
    auprc = average_precision_score(y_true, y_scores)
    
    fpr, tpr, thresholds = roc_curve(y_true, y_scores)
    
    # FPR @ 95% TPR
    idx = np.where(tpr >= 0.95)[0][0]
    fpr95 = fpr[idx]
    
    return {
        'auroc': float(auroc),
        'auprc': float(auprc),
        'fpr95': float(fpr95)
    }

def compute_calibration_metrics(
    map_est: torch.Tensor, 
    uncertainty: torch.Tensor, 
    gt: torch.Tensor
) -> Dict[str, float]:
    """
    Compute calibration metrics: Coverage, Sharpness, and standardized residual z.
    Uncertainty here is IQR, which corresponds to roughly 1.349 * sigma for a Gaussian.
    We'll use it to define the CI bounds.
    """
    map_est_np = map_est.numpy()
    uncertainty_np = uncertainty.numpy()
    gt_np = gt.numpy()
    
    # We use uncertainty as width of the interval for a simple coverage heuristic
    # (Assuming uncertainty is the width of a ~50% or 95% CI depending on config, 
    # uGUIDE returns IQR by default, so it's a 50% CI. We'll measure 50% coverage).
    lower_bound = map_est_np - (uncertainty_np / 2.0)
    upper_bound = map_est_np + (uncertainty_np / 2.0)
    
    coverage = np.mean((gt_np >= lower_bound) & (gt_np <= upper_bound), axis=0)
    sharpness = np.mean(uncertainty_np, axis=0)
    
    # Standardized residual: z = (MAP - GT) / (IQR / 1.349)
    # Adding epsilon to prevent division by zero
    sigma_approx = (uncertainty_np / 1.349) + 1e-8
    z_score = (map_est_np - gt_np) / sigma_approx
    
    return {
        'coverage_50_mean': float(np.mean(coverage)),
        'sharpness_mean': float(np.mean(sharpness)),
        'z_score_magnitude': np.abs(z_score) # Used for coupling analysis
    }

def compute_coupling(ood_scores: np.ndarray, z_score_magnitude: np.ndarray) -> float:
    """
    Compute the Spearman correlation between OOD scores and calibration failure.
    Calibration failure is measured by the magnitude of the standardized residual.
    We average the z-score magnitude over the parameter dimensions.
    """
    # Average z_score magnitude across parameters for an overall unreliability scalar
    mean_z_mag = np.mean(z_score_magnitude, axis=1)
    
    rho, _ = spearmanr(ood_scores, mean_z_mag)
    return float(rho)
