
# BatchMap - A parallel implementation of the OneMap R package for fast computation of linkage maps in outcrossing species

BatchMap is a fork of the [OneMap](https://github.com/augusto-garcia/onemap) software package for the construction of linkage maps.
Specifically it aims to provide a framework that enables creating linkage maps from dense marker data (N>10,000).

Most of the non-core functionality of OneMap has been stripped and only a few functions are left at user-level. Further, this package is typically expected to be run on a headless server due to memory requirements, so all graphic functionality has also been removed.

## Why BatchMap

BatchMap employs a divide and conquer algorithm that manages to speed up the calculation of high density datasets and - additionally - scales well with higher marker numbers. It further features a routine inspired by OneMap's `ripple.seq` function that can adaptively reorder windows of markers to improve the order of sequences.

## Installation
### Install from GitHub

```R
install.packages("devtools")
library(devtools)
install_github("bschiffthaler/BatchMap")
```

## Vignette - tutorial

After installation, please type `vignette("BatchMap")` to open the long format documentation.

## Citation

This work depends whole-heartedly on the work done by all original and current authors of OneMap. We therefore ask you to cite both OneMap and BatchMap if you find BatchMap useful in your research.
