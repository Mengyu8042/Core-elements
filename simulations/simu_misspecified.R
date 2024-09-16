rm(list = ls())
gc()
## Import R Packages ##
list.of.packages <- c("mvtnorm", "MASS", "corpcor", "pracma", "Rfast",
                      "Rcpp", "Matrix", "ggplot2", "gridExtra", "cowplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)
library(mvtnorm)
library(MASS)
library(corpcor)
library(pracma)
library(Rfast)
library(Rcpp)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(cowplot)

## Import the source file ##
setwd(getwd())
source("../source/utils.R")
sourceCpp("../source/funcs.cpp")
sourceCpp("../source/oss/myoss.cpp")

## Initialization ##
set.seed(1000)
is_predict <- FALSE

file_name <- ifelse(is_predict, "pmse_misspecified", "mse_misspecified")

method_names <- c("FULL", "UNIF", "DOUBLY", "BLEV", "SLEV", 
                  "IBOSS", "OSS", "DOPT", "CORE")
num_method <- length(method_names)
method_colors <- c("#1a1a1a", "#a6cee3", "#1f78b4", "#cab2d6", "#6a3d9a",     
                   "#fa9fb5", "#ae017e", "#bf812d", "red")


fig_all <- list()
fig_num <- 0

## Beginning ##
for (noise_type in c("M1", "M2", "M3")) {
  for (spar_type in c("R1", "R2", "R3", "R4", "R5")) {
    fig_num <- fig_num + 1
    case <- c("D1", spar_type)
    print(c(fig_num, case))

    ## Set parameters ##
    N <- 10000  # sample size n
    p <- 20  # dimension p
    sub_meta <- (1:5) * 2 * p  # subsample parameter r
    nloop <- 100  # number of replicates
    SNR <- 4  # signal-to-noise ratio
    beta_true <- rep(1, p)  # true coefficients
    
    mse_meta <- array(0, dim = c(nloop, num_method, length(sub_meta)))
    mse_temp <- matrix(0, num_method, length(sub_meta))

    ##########################################################################
    for(i in 1:nloop){
      
      ## Set the covariance matrix of X ##
      Sigma <- diag(p)
      for (ii in 1:p) {
        for (jj in 1:p) {
          Sigma[ii, jj] <- 0.6^(abs(ii - jj))
        }
      }
      set.seed(100 + 123 * i)

      ## Generate synthetic data ##
      if (!is_predict) {
        data <- dataGenerate(N, p, Sigma, beta_true, SNR, case, 
                             is_corrupt = FALSE, n_bad = 0, noise_type = noise_type)
        X <- data$X
        Y <- data$Y
        
      } else {
        N_train <- floor(0.7 * N)
        N_test <- N - N_train
        
        data <- dataGenerate(N_train, p, Sigma, beta_true, SNR, case, 
                             is_corrupt = FALSE, n_bad = 0, noise_type = noise_type)
        X <- data$X
        Y <- data$Y
        
        data <- dataGenerate(N_test, p, Sigma, beta_true, SNR, case, 
                             is_corrupt = FALSE, n_bad = 0, noise_type = noise_type)
        X_test <- data$X
        Y_test <- data$Y
      }
      rm(data)

      ## Full sample results ##
      res <- method_FULL(X, Y)
      full_temp <- res$beta_est
      prob_blev <- res$prob_blev
      prob_slev <- res$prob_slev
      ord_dopt <- res$ord_dopt
      
      for(j in 1:length(sub_meta)){
        set.seed(100 + 421 * j + 123 * i)
        sub <- sub_meta[j]
        
        unif_temp <- method_UNIF(X, Y, sub)
        doubly_temp <- method_DOUBLY(X, Y, sub)
        blev_temp <- method_BLEV(X, Y, sub, prob_blev)
        slev_temp <- method_SLEV(X, Y, sub, prob_slev)
        iboss_temp <- method_IBOSS(X, Y, sub)
        oss_temp <- method_OSS(X, Y, sub)
        dopt_temp <- method_DOPT(X, Y, sub, ord_dopt)
        core_temp <- method_CORE(X, Y, sub)
        
        
        ######################################################
        coef <- cbind(full_temp, unif_temp, doubly_temp, blev_temp, slev_temp, 
                      iboss_temp, oss_temp, dopt_temp, core_temp)
        coef[is.na(coef)] <- 0
        if (!is_predict) {
          coef_mse <- colMeans((coef - beta_true)^2)/mean(beta_true^2)
        } else {
          coef_mse <- colMeans((X_test %*% coef - Y_test)^2)/mean(Y_test^2)
        }
        mse_temp[, j] <- coef_mse
      }
      
      mse_meta[i, , ] <- mse_temp
    }
    
    ## Plot the results ##
    log_mse_list <- vector("list", num_method)
    for (ii in 1:num_method) {
      log_mse_list[[ii]] <- log10(mse_meta[, ii, ])
    }
    
    mse_mat <- data.frame(
      mse = unlist(lapply(log_mse_list, function(x) apply(x, 2, mean))),
      sd = unlist(lapply(log_mse_list, function(x) apply(x, 2, sd))),
      Method = factor(rep(method_names, each = length(sub_meta))),
      sub = rep(sub_meta/p, num_method)
    )
    mse_mat$Method <- factor(mse_mat$Method, levels = method_names)
    
    fig_title <- paste(c(noise_type, spar_type), collapse = ", ")
    p23 <- plot_mse(mse_mat, sub_meta, is_predict, method_colors, fig_title)
    print(p23)
    fig_all[[fig_num]] <- p23
    
    # save(mse_mat, sub_meta, is_predict, method_colors, fig_title, 
    #      file = paste0(file_name, "_", noise_type, "_", spar_type, ".RData"))
  }
}

legend <- get_legend(p23)
fig_all <- lapply(fig_all, function(x) x + theme(legend.position = "none"))
combined_plots <- plot_grid(plotlist = fig_all, nrow = 3, ncol = 5)
final_plot <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(10, 1))

# save(final_plot, file = paste0(file_name, ".RData"))

jpeg(paste0(file_name, ".jpg"), width = 17.5, height = 11.5, units = "in", res = 400)
print(final_plot)
dev.off()
