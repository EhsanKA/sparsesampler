# SParseSampler (SPS)

[![PyPI version](https://badge.fury.io/py/sparsesampler.svg)](https://badge.fury.io/py/sparsesampler)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

SParseSampler (SPS) is a Python package for efficient subsampling of large-scale single-cell RNA-seq and flow cytometry datasets while preserving rare cell populations. The method employs an unsupervised approach that maintains dataset structure and rare cell types without requiring explicit labels.

## Key Features

### Core Method
- PCA-based dimensionality reduction with automatic parameter selection via an EVR-based heuristic
- Variance-weighted binning in the reduced dimensional space
- Iterative cell selection prioritizing sparsest bins
- Preserves rare populations without requiring cell type labels

### Performance Benefits
- Computational efficiency comparable to random sampling
- Superior rare cell retention compared to existing methods (Hopper, Atomic Sketch)
- Performance comparable to scSampler with improved speed
- Successfully tested on datasets up to 34 million cells
- Validated at multiple rarity levels (1%, 0.5%, 0.1%)

## Installation

```bash
pip install sparsesampler
```

## Technical Details

### Parameters

SParseSampler uses an **EVR-based heuristic** for automatic parameter selection. Both the number of principal components (`p`) and the Bin Resolution Factor (`K`) are derived from the explained variance ratio (EVR) of the principal components using the EVR index (`feature_index`).

- **EVR Index (`feature_index`)**
  - Default: 12
  - Controls both the number of principal components and the bin resolution
  - The number of principal components `p` is set to `feature_index + 1`
  - The Bin Resolution Factor `K` is computed as `K = 2 / EVR_i`, where `EVR_i` is the explained variance ratio of the `i`-th principal component
  - **Recommended ranges:**
    - Flow cytometry data: EVR indices 7–20
    - scRNA-seq data: EVR indices 12–25
  - The default value of 12 lies within both optimal ranges

- **Sample Size (`size`)**
  - Number of cells to subsample from the dataset

- **Seed (`seed`)**
  - Random seed for reproducibility

### Supported Data Types
- Single-cell RNA sequencing data
- Flow cytometry data

## Validation

### Benchmarking
- Comprehensive comparison against state-of-the-art methods
- Validated on large-scale datasets:
  - MCC dataset (scRNA-seq): 3.2M cells, 3,065 genes
  - LCMV dataset (flow cytometry): 34M cells, 31 features
- Consistent performance across varying dataset sizes and rarity levels (1%, 0.5%, 0.1%)
- Downstream validation: Random Forest classifiers trained on SPS-subsampled data achieve substantially higher F1 scores for rare cell types compared to random subsampling

### PCA Runtime
- Flow cytometry (LCMV, 34M cells, 31 features): ~11 seconds
- scRNA-seq (MCC, 3M cells, 3,065 genes): ~5 minutes

## Usage

```python
import sparsesampler.sampling as sps
import numpy as np

# Load your data (n_samples × n_features)
# Example 1: From NumPy array
X = np.load('your_data.npy')

# Example 2: From CSV file
import pandas as pd
X = pd.read_csv('your_data.csv').values

# Example 3: scRNA-seq data (AnnData format)
import scanpy as sc
adata = sc.read_h5ad('your_data.h5ad')
X = adata.X  # Use .toarray() if sparse matrix

# Run SParseSampler with default parameters (EVR index = 12)
indices, elapsed_time = sps.sample(X=X, size=100000)

# Run with custom EVR index (e.g., for flow cytometry data)
indices, elapsed_time = sps.sample(X=X, size=50000, feature_index=8)

# Get subsampled data
X_sampled = X[indices]
```

### Preprocessing Recommendations

We recommend applying standard quality control filtering prior to SPS, including:
- Removal of cells with abnormally high/low UMI counts
- Filtering cells with high mitochondrial gene percentages
- Doublet detection and removal (e.g., using Scrublet or DoubletFinder)

## Visualization

The following animation shows how points are selected from a 2D toy dataset using PCA binning. Points are selected category by category (cells with 1 point, 2 points, etc.), and the process is visualized step by step:

- All points start as skyblue.
- When a category is considered, the cells are highlighted in yellow and the points in those cells are shown in gray for visibility.
- Selected points turn red and remain red in all subsequent frames.
- The process continues until the target number of points is reached.

![Sampling Process Animation](./docs/sampling_process.gif)

To generate the animation yourself, run:

```bash
python docs/generate_visualization.py
```

## Citation

If you use SParseSampler in your research, please cite:

```bibtex
# Add citation when available
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
