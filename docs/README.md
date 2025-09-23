# Documentation and Visualization

This directory contains files used for generating documentation and visualizations for the SParseSampler project. **These files are not part of the main package** and are intended solely for creating materials for the GitHub repository and documentation.

## Contents

- `generate_visualization.py`: Script to generate the animated sampling process demonstration
- `sampling_process.gif`: Animated visualization showing how SParseSampler works
- `sampling_process_frames/`: Individual frames from the animation

## Usage

To regenerate the visualization:

```bash
python docs/generate_visualization.py
```

This will create:
- `docs/sampling_process.gif`: The main animation
- `docs/sampling_process_frames/`: Directory with individual frame images

## Purpose

The visualization demonstrates the step-by-step sparse sampling process using toy data and PCA binning to help users understand how the algorithm works:

1. Shows original 2D toy data
2. Displays PCA projection with binning grid
3. Iteratively highlights and selects cells with the fewest points
4. Shows the final sampled subset

This is purely for educational and documentation purposes and should not be used in production code.