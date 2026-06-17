import torch
import numpy as np

class MahalanobisDetector:
    """
    Family 1 Detector: Summary-space OOD detector using Mahalanobis distance.
    Fits the ID distribution of summary embeddings and scores test embeddings.
    """
    def __init__(self):
        self.mean = None
        self.inv_cov = None

    def fit(self, id_embeddings: torch.Tensor):
        """
        Fit the detector on ID summary embeddings.
        
        Parameters
        ----------
        id_embeddings : torch.Tensor
            Embeddings of ID samples (shape: [N, D]).
        """
        id_embeddings_np = id_embeddings.cpu().numpy()
        self.mean = np.mean(id_embeddings_np, axis=0)
        cov = np.cov(id_embeddings_np, rowvar=False)
        # Add small regularization to diagonal to prevent singular matrix
        cov += np.eye(cov.shape[0]) * 1e-6
        self.inv_cov = np.linalg.inv(cov)

    def score(self, test_embeddings: torch.Tensor) -> torch.Tensor:
        """
        Compute the OOD score (Mahalanobis distance) for test embeddings.
        Higher score means more OOD.
        
        Parameters
        ----------
        test_embeddings : torch.Tensor
            Embeddings of test samples (shape: [N, D]).
            
        Returns
        -------
        torch.Tensor
            OOD scores (shape: [N]).
        """
        assert self.mean is not None, "Detector must be fitted before scoring."
        
        test_embeddings_np = test_embeddings.cpu().numpy()
        diff = test_embeddings_np - self.mean
        # Mahalanobis distance: sqrt( (x - u)^T * Cov^-1 * (x - u) )
        dist_sq = np.sum(np.dot(diff, self.inv_cov) * diff, axis=1)
        dist = np.sqrt(dist_sq)
        
        return torch.tensor(dist, dtype=torch.float32)
