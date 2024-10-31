import random
import time
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sparseFlow.preprocessing import perform_pca_binning, adjust_feature_importances, accumulate_indices_until_threshold


def generate_sparse_sample(adata, size, seed=1234):
    scaler = StandardScaler()
    data_standardized = scaler.fit_transform(adata.X)
    X = data_standardized

    random.seed(seed)
    print(f'********* #Start# *********')
    start_time = time.time()

    n_components = min(adata.shape[1], 100)
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
