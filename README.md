# Core-Elements for Classical Linear Regression
This repository includes the implementation of our work **"Core-Elements for Classical Linear Regression"**.


## Introduction
A brief introduction about the folders and files:
* `data/`: open-source datasets used in Section 6;
* `source/`: the source file for sparse matrix operations;
* `simu_estimation.R`, `simu_prediction.R`: simulation scripts for estimation and prediction, respectively;
* `example1_estimation.R`, `example1_prediction.R`: data analysis scripts for Example 1 in Section 6;
* `example2_estimation.R`, `example2_prediction.R`: data analysis scripts for Example 2 in Section 6;
* `demo.Rmd`, `demo.html`: a simulation example under the case of **D1, R1**.


## Reproducibility
For simulation study in Section 5,
* you can run `simu_estimation.R` with `is_corrupt <- FALSE` (or `is_corrupt <- TRUE`) to reproduce the results in Fig. 2 (or Fig. 5);
* you can run `simu_prediction.R` with `is_corrupt <- FALSE` (or `is_corrupt <- TRUE`) to reproduce the results in Fig. 3 (or Fig. 6).

For data analysis in Section 6,
* you can run `example1_estimation.R` and `example1_prediction.R` to reproduce the results in Fig. 7;
* you can run `example2_estimation.R` and `example2_prediction.R` to reproduce the results in Fig. 8.