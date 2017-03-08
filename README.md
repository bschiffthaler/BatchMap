
# BatchMap - A parallel implementation of the OneMap R package for fast computation of F1 linkage maps in outcrossing species

BatchMap is a fork of the [OneMap](https://github.com/augusto-garcia/onemap) software package for the construction of linkage maps.
Specifically it aims to provide a framework that enables creating linkage maps from dense marker data (N>10,000).

Most of the non-core functionality of OneMap has been stripped and only a few functions are left at user-level. Further, this package is typically expected to be run on a headless server due to memory requirements, so all graphic functionality has also been removed.

## Why use BatchMap

BatchMap employs a divide and conquer algorithm that manages to speed up the calculation of high density datasets and - additionally - scales well with higher marker numbers. It further features a routine inspired by OneMap's `ripple.seq` function that can adaptively reorder windows of markers to improve the order of sequences.

## When not to use BatchMap

BatchMap is created specifically for F1 outcrossing populations. If your data is not that (e.g. backcross, F2, RIL), use [OneMap](https://github.com/augusto-garcia/onemap).

## Memory requirements

The twopoint table of recombination fractions and likelihoods will take `4 * 8 * N * N bytes` in memory, where `N` is your marker number. If you have anything less than that, you need to either reduce the number of input markers, or rent a cloud server (Amazon EC2 provides high memory machines). 

## Installation
### Install from GitHub

```R
install.packages("devtools")
library(devtools)
install_github("bschiffthaler/BatchMap")
```

## Note for Windows users

Windows is currently not able to parallelize R code using the `parallel` package. It is therefore **highly** recommended to use the 

## Vignette - tutorial

After installation, you can read the user tutorial ('vignette') [here]().

## Citation

This software depends whole-heartedly on the work done by all original and current authors of OneMap. We therefore ask you to cite both OneMap and BatchMap if you find BatchMap useful in your research.
