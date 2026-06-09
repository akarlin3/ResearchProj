def run_tier3(data_path: str, model_config_path: str):
    """
    Tier 3 Experiment: Application to glioma dataset.
    This tier tests performance on a different clinical cohort and evaluates downstream clinical impact.
    
    Parameters
    ----------
    data_path : str
        Path to the glioma dataset.
    model_config_path : str
        Path to the trained model configuration.
    """
    raise NotImplementedError("Tier 3 is out of scope for the current build.")

if __name__ == "__main__":
    run_tier3("./data/glioma", "./results/config.pkl")
