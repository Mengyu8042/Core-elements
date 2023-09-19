# Core-Elements for Large-Scale Least Squares Estimation
This repository includes the implementation of our work **"Core-Elements for Large-Scale Least Squares Estimation"**.


## Introduction
A brief introduction about the folders and files:
* `data/`: the real-world dataset used in Section 6;
* `source/`: the source file for sparse matrix operations;
* `simu_*.R`: simulation scripts for Section 5 and Supplementary material;
* `realdata_*.R`: real data analysis scripts for Section 6.


## Reproducibility
For simulation studies in Section 5,
* you can run `simu_estimation.R` with `is_corrupt <- FALSE` (resp., `is_corrupt <- TRUE`) to reproduce Fig. 2 (resp., Fig. 5);
* you can run `simu_prediction.R` with `is_corrupt <- FALSE` (resp., `is_corrupt <- TRUE`) to reproduce Fig. 3 (resp., Fig. 6).
* you can run `simu_error_bound.R` to reproduce Fig. 4.

For real data analysis in Section 6,
* you can run `realdata_visualization.R` to reproduce Fig. 7;
* you can run `realdata_estimation.R` and `realdata_prediction.R` to reproduce Fig. 8.

For additional numerical results in Section 2 of the supplementary material,
* you can run `simu_estimation_misspecified.R` to reproduce Fig. S1;
* you can run `simu_prediction_misspecified.R` to reproduce Fig. S2.
