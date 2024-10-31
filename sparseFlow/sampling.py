import random
import time
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sparseFlow.preprocessing import perform_pca_binning, adjust_feature_importances, accumulate_indices_until_threshold


def sample(adata=None, X=None, size=50000, seed=1234):
    """
    Perform PCA and binning to sample cells based on the PCA space.
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    X: np.ndarray
        Data matrix.
    size: int
        Number of cells to sample.
    seed: int
        Random seed.
    Returns
    -------
    samples: list
        List of indices of sampled cells.
    elapsed_time: float
        Elapsed time in seconds.
    Raises
    ------
    ValueError
        Either adata or X must be provided. If both are provided, adata will be used.
    """

    scaler = StandardScaler()
    if adata is not None:
        data_standardized = scaler.fit_transform(adata.X)
    elif X is not None:
        data_standardized = scaler.fit_transform(X)
    else:
        raise ValueError("Either adata or X must be provided.")

    X = data_standardized

    random.seed(seed)
    print(f'********* #Start# *********')
    start_time = time.time()

    n_components = min(X.shape[1], 100)
    pca = PCA(n_components=n_components)
    pca.fit(X)
    X_pca = pca.transform(X)
    df = pd.DataFrame(X_pca[:, :n_components], columns=[f'PC{i + 1}' for i in range(n_components)])

    feature_importances = adjust_feature_importances(pca)
    perform_pca_binning(df, feature_importances)

    threshold = size
    samples = accumulate_indices_until_threshold(df, threshold, seed=seed)

    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")
    
    return samples, elapsed_time
