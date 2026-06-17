import torch
import numpy as np
import json
from pathlib import Path

from sibyl.data.synthetic import generate_id_dataset, PRIOR_BOUNDS
from sibyl.data.shift import generate_shifted_datasets
from sibyl.estimator.wrapper import create_and_train_estimator, get_posterior, get_embedding, test_embedding_matches_uguide
from sibyl.detectors.family1 import MahalanobisDetector
from sibyl.detectors.family2 import ResidualConformalDetector
from sibyl.metrics.eval import compute_detection_metrics, compute_calibration_metrics, compute_coupling
from sibyl.forward_model.ivim import ivim_biexponential, DEFAULT_B_SCHEME

def run_tier1(folderpath: str = './results'):
    folderpath_p = Path(folderpath)
    folderpath_p.mkdir(parents=True, exist_ok=True)
    
    # 1. Generate Training Data
    print("Generating synthetic datasets...")
    N_TRAIN = 10000
    N_CALIB = 2000
    N_TEST = 1000
    
    theta_train, x_train = generate_id_dataset(N_TRAIN, snr=50.0, seed=42)
    theta_calib, x_calib = generate_id_dataset(N_CALIB, snr=50.0, seed=43)
    theta_test_id, x_test_id = generate_id_dataset(N_TEST, snr=50.0, seed=44)
    
    # Generate clean signal for shift injection on test set
    f, D, Dstar = theta_test_id[:, 0], theta_test_id[:, 1], theta_test_id[:, 2]
    clean_test_signal = ivim_biexponential(DEFAULT_B_SCHEME, f, D, Dstar)
    
    # 2. Train Estimator
    config = create_and_train_estimator(theta_train, x_train, folderpath=folderpath, epochs=100, seed=42)
    
    # 3. Smoke Test
    print("Running Smoke Test...")
    test_embedding_matches_uguide(x_train[:5], config)
    
    # 4. Prepare Detectors (Calibration)
    print("Fitting detectors...")
    # Family 1 (Summary-space)
    calib_embeddings = get_embedding(x_calib, config)
    det1 = MahalanobisDetector()
    det1.fit(calib_embeddings)
    
    # Family 2 (Signal-space / Residuals)
    # We need MAP estimate for the calib set to get the predicted signal
    map_calib, _ = get_posterior(x_calib, config)
    f_c, D_c, Dstar_c = map_calib[:, 0], map_calib[:, 1], map_calib[:, 2]
    pred_signal_calib = ivim_biexponential(DEFAULT_B_SCHEME, f_c, D_c, Dstar_c)
    calib_residuals = x_calib - pred_signal_calib
    
    det2 = ResidualConformalDetector()
    det2.fit(calib_residuals)
    
    # 5. Evaluate ID Test Set
    id_embeddings = get_embedding(x_test_id, config)
    map_id, unc_id = get_posterior(x_test_id, config)
    pred_signal_id = ivim_biexponential(DEFAULT_B_SCHEME, map_id[:, 0], map_id[:, 1], map_id[:, 2])
    id_residuals = x_test_id - pred_signal_id
    
    scores_id_det1 = det1.score(id_embeddings).numpy()
    scores_id_det2 = det2.score(id_residuals).numpy()
    
    # Control: Degeneracy-FPR
    # Define intrinsically degenerate samples: e.g., extremely low f (f < 0.05), where D* is poorly identified.
    degenerate_mask = theta_test_id[:, 0] < 0.05
    n_degen = degenerate_mask.sum().item()
    
    # Thresholds for FPR based on 95th percentile of ALL ID scores (5% baseline FPR)
    thresh_det1 = np.percentile(scores_id_det1, 95)
    thresh_det2 = np.percentile(scores_id_det2, 95)
    
    fpr_degen_det1 = np.mean(scores_id_det1[degenerate_mask.numpy()] > thresh_det1)
    fpr_degen_det2 = np.mean(scores_id_det2[degenerate_mask.numpy()] > thresh_det2)
    
    print(f"\n--- Degeneracy FPR Control (Baseline is 0.05) ---")
    print(f"Number of degenerate ID samples (f < 0.05): {n_degen}")
    print(f"Family 1 (Summary) FPR on degenerate: {fpr_degen_det1:.4f}")
    print(f"Family 2 (Residual) FPR on degenerate: {fpr_degen_det2:.4f}")
    
    # 6. Evaluate Shifted Datasets
    shifted_datasets = generate_shifted_datasets(theta_test_id, clean_test_signal, seed=45)
    
    results = {
        'degeneracy_control': {
            'det1_fpr': float(fpr_degen_det1),
            'det2_fpr': float(fpr_degen_det2)
        },
        'shifts': {}
    }
    
    print("\n--- Shift Evaluation ---")
    for shift_name, (theta_ood, x_ood) in shifted_datasets.items():
        # Extact OOD variables
        ood_embeddings = get_embedding(x_ood, config)
        map_ood, unc_ood = get_posterior(x_ood, config)
        pred_signal_ood = ivim_biexponential(DEFAULT_B_SCHEME, map_ood[:, 0], map_ood[:, 1], map_ood[:, 2])
        ood_residuals = x_ood - pred_signal_ood
        
        scores_ood_det1 = det1.score(ood_embeddings).numpy()
        scores_ood_det2 = det2.score(ood_residuals).numpy()
        
        # Detection metrics
        det1_metrics = compute_detection_metrics(scores_id_det1, scores_ood_det1)
        det2_metrics = compute_detection_metrics(scores_id_det2, scores_ood_det2)
        
        # Calibration metrics on the OOD set
        calib_metrics = compute_calibration_metrics(map_ood, unc_ood, theta_ood)
        
        # Coupling metrics
        # Combine ID and OOD for coupling so we see the correlation over the full distribution of trustworthiness
        all_scores_det1 = np.concatenate([scores_id_det1, scores_ood_det1])
        all_scores_det2 = np.concatenate([scores_id_det2, scores_ood_det2])
        
        id_calib = compute_calibration_metrics(map_id, unc_id, theta_test_id)
        all_z_mag = np.concatenate([id_calib['z_score_magnitude'], calib_metrics['z_score_magnitude']])
        
        coupling_det1 = compute_coupling(all_scores_det1, all_z_mag)
        coupling_det2 = compute_coupling(all_scores_det2, all_z_mag)
        
        results['shifts'][shift_name] = {
            'det1_metrics': det1_metrics,
            'det2_metrics': det2_metrics,
            'calib_metrics': {k: v for k, v in calib_metrics.items() if k != 'z_score_magnitude'},
            'coupling_rho_det1': coupling_det1,
            'coupling_rho_det2': coupling_det2
        }
        
        print(f"\n[{shift_name}]")
        print(f"  Det1 AUROC: {det1_metrics['auroc']:.3f} | Det2 AUROC: {det2_metrics['auroc']:.3f}")
        print(f"  Det1 Coupling (Rho): {coupling_det1:.3f} | Det2 Coupling (Rho): {coupling_det2:.3f}")
        print(f"  Coverage (50% CI): {calib_metrics['coverage_50_mean']:.3f}")

    # Save results
    out_file = folderpath_p / 'tier1_results.json'
    with open(out_file, 'w') as f:
        json.append = False # Just ensuring standard
        json.dump(results, f, indent=4)
        
    print(f"\nSaved results to {out_file}")

if __name__ == "__main__":
    run_tier1()
