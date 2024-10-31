import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import time
import scipy.sparse


def pca_bin_sample_(df, feature_importances, seed=12345):
    np.random.seed(seed)

    # Ensure num_pcs does not exceed the number of columns in df
    num_pcs = min(feature_importances.shape[0], df.shape[1])

    # Function to create bins and digitize
    def create_bins_and_digitize(data, n_bins):
        edges = np.linspace(data.min(), data.max(), n_bins + 1)
        bins = np.digitize(data, edges)
        return bins

    def compute_sample_bins(df, bin_sizes):
        bins = [create_bins_and_digitize(df.iloc[:, i], bin_sizes[i]) for i in range(num_pcs)]

        # Combine bins to form grid cells
        df['grid_cell'] = list(zip(*bins))
        
        return 

    compute_sample_bins(df, feature_importances)
    return

    
def set_min_to_two(pca):
    if len(pca.explained_variance_ratio_)< 40:
        out = np.ceil(pca.explained_variance_ratio_*100).astype(int)

    else:
        mul = 2.0/ pca.explained_variance_ratio_[19]
        out = np.ceil(pca.explained_variance_ratio_*mul).astype(int)
        
    return out[out>2]


import random

def find_threshold_index(sorted_grid_cells, threshold):
    cumulative = 0
    for index, frequency in sorted_grid_cells.value_counts().sort_index().items():
        cumulative += index * frequency
        if cumulative >= threshold:
            return index
    return None

def accumulate_indices_until_threshold(df, threshold, seed=1234):
    random.seed(seed)
    # Count the occurrences of each grid_cell
    grid_cell_counts = df['grid_cell'].value_counts()

    # Sort the grid_cells by count in ascending order
    sorted_grid_cells = grid_cell_counts.sort_values()

    # Find the threshold index
    threshold_index = find_threshold_index(sorted_grid_cells, threshold)
    print(f'threshold_index is : {threshold_index}')
    
    # Group the DataFrame by 'grid_cell'
    grouped_df = df.groupby('grid_cell')

    accumulated_indices = []
    accumulated_count = 0
    all_remainings_indices = []

    # Iterate over sorted grid_cells and accumulate indices
    for grid_cell in sorted_grid_cells.index:
        group_indices = grouped_df.get_group(grid_cell).index.tolist()
        if len(group_indices) < threshold_index:
            accumulated_indices.extend(group_indices)
            accumulated_count += len(group_indices)
        elif len(group_indices) == threshold_index:
            all_remainings_indices.extend(group_indices)
        else:
            break


    # Calculate how many more indices we need to reach the threshold
    remaining_count = threshold - accumulated_count
    print(f'remaining is : {remaining_count}')

    # Randomly select the remaining indices from the current group
    accumulated_indices.extend(random.sample(all_remainings_indices, remaining_count))
    
    return accumulated_indices

def generate_cubic(adata, size, seed = 1234):
    
    scaler = StandardScaler()
    data_standardized = scaler.fit_transform(adata.X)
    X = data_standardized


    random.seed(seed)
    print(f'********* #Start# *********')
    start_time = time.time()

    N_components= min(adata.shape[1], 100)

    pca = PCA(n_components=N_components)
    pca.fit(X)
    X_pca = pca.transform(X)
    df = pd.DataFrame(X_pca[:, :N_components], columns=[f'PC{i+1}' for i in range(N_components)])


    feature_importances = set_min_to_two(pca)

    pca_bin_sample_(df, feature_importances)

    threshold = size  # Set your desired threshold
    samples = accumulate_indices_until_threshold(df, threshold, seed=seed)

    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

        
    return samples, elapsed_time

