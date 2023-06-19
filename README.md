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
* you can run `example1_estimation.R` and `example1_prediction.R` to reproduce the results in Fig. 8;
* you can run `example2_estimation.R` and `example2_prediction.R` to reproduce the results in Fig. 10.

For data analysis in Section S2.2,
* you can run `example1_estimation.R` and `example1_prediction.R` (resp. `example2_estimation.R` and `example2_prediction.R`) with `boot_type <- "boot_pair"` to reproduce the results in Fig. 3(a) (resp. Fig. 3(b));
* you can run `example2_estimation.R` and `example2_prediction.R` with `simu_y <- TRUE` and `boot_type <- "boot_pair"` (resp. `boot_type <- "boot_res"`) to reproduce the results in Fig. 5(a) (resp. Fig. 5(b)).