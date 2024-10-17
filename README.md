# Core-Elements for Large-Scale Least Squares Estimation
This repository includes the implementation of our work **"Core-Elements for Large-Scale Least Squares Estimation"** [https://arxiv.org/abs/2206.10240].


If you use this toolbox in your research and find it useful, please cite:
```
@article{li2024core,
  title={Core-elements for large-scale least squares estimation},
  author={Li, Mengyu and Yu, Jun and Li, Tao and Meng, Cheng},
  journal={Statistics and Computing},
  volume={34},
  number={6},
  pages={1--16},
  year={2024},
  publisher={Springer}
}
```


## Introduction
A brief introduction about the folders and files:

* `simulations/`: simulation scripts;
  - `simu_main.R`: reproduces the numerical results in Section 5;
  - `simu_misspecified.R`: reproduces the numerical results in the Appendix.
* `source/`: source files;
  - `oss/`: implementation code for the OSS method provided by its authors;
  - `funcs.cpp`: useful functions for sparse matrix computation;
  - `utils.R`: useful functions for data generation, method implementation, plotting, and more.


## Reproducibility

### Simulation Studies (Section 5)
* To reproduce the simulation results, you can run `simu_main.R`:
  - Set `is_corrupt <- FALSE` and `is_predict <- FALSE` to generate Figure 2.
  - Set `is_corrupt <- FALSE` and `is_predict <- TRUE` to generate Figure 3.
  - Set `is_corrupt <- TRUE` and `is_predict <- FALSE` to generate Figure 5.
  - Set `is_corrupt <- TRUE` and `is_predict <- TRUE` to generate Figure 6.


### Additional Numerical Results (Appendix Section B)
* To reproduce additional results:
  - Run `simu_misspecified.R` with `is_predict <- FALSE` to generate Figure 1.
  - Run `simu_misspecified.R` with `is_predict <- TRUE` to generate Figure 2.