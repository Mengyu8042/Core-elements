rm(list = ls())
gc()
## Import R Packages ##
list.of.packages <- c("mvtnorm", "MASS", "corpcor", "pracma", "Rfast",
                      "Rcpp", "Matrix", "ggplot2", "gridExtra")
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

## Import the source file ##
setwd(getwd())
sourceCpp("source/funcs.cpp")

## Initialization ##
set.seed(1000)
fig_all <- list()
fig_num <- 0
is_corrupt <- FALSE  # add outliers or not
if (is_corrupt == TRUE) {
  file_name <- "mse_corrupt.pdf"
} else {
  file_name <- "mse_uncorrupt.pdf"
}

## Define some functions ##
yGenerate <- function(x, beta, SNR) {
  # generate the response vector
  y_raw <- as.vector(x %*% beta)
  N <- length(y_raw)
  y <- y_raw + rnorm(N, 0, sd(y_raw)/sqrt(SNR))
  return(y)
}

xSparse <- function(x, sparsity) {
  # sparsify the predictor matrix
  N <- length(x[, 1])
  p <- length(x[1, ])
  xx <- rep(0, N * p)
  x_vec <- as.vector(x)
  id <- sample(1:(N * p), (1 - sparsity) * N * p, replace = FALSE)  # non-zeros
  xx[id] <- x_vec[id]
  xx <- matrix(xx, N, p)
  return(xx)
}

dataCorrupt <- function(n_bad, x, y, beta, SNR) {
  # corrupt the iid data
  N <- length(x[, 1])
  p <- length(x[1, ])
  
  m <- ceil(n_bad/4)
  x_bad1 <- -10 + matrix(rnorm(m * p), m, p)
  y_bad1 <- 1000 + 10 * rnorm(m)
  x_bad2 <- 10 + matrix(rnorm(m * p), m, p)
  y_bad2 <- -500 + 10 * rnorm(m)
  x_bad3 <- matrix(runif(m * p), m, p)
  y_bad3 <- rbinom(m, 1, 0.5)
  x_bad4 <- mvrnorm(n_bad - 3 * m, rep(0, p), diag(p))
  y_raw <- as.vector(x_bad4 %*% beta)
  y_bad4 <- y_raw + rt(n_bad - 3 * m, df = 2)
  x_corrupt <- rbind(x[1:(N - n_bad), ], x_bad1, x_bad2, x_bad3, x_bad4)
  y_corrupt <- c(y[1:(N - n_bad)], y_bad1, y_bad2, y_bad3, y_bad4)
  
  return(list(x_corrupt, y_corrupt))
}


## Beginning ##
for (case_type in c("D1", "D2", "D3")) {
  for (spar_type in c("R1", "R2", "R3", "R4", "R5")) {
    fig_num <- fig_num + 1
    case <- c(case_type, spar_type)
    print(c(fig_num, case))
    
    ## Set parameters ##
    N <- 10000  # sample size n
    if (is_corrupt == TRUE) {
      p <- 20  # dimension p
      sub_meta <- (1:5) * 5 * p  # subsample parameter r
    } else {
      p <- 100
      sub_meta <- (1:5) * 2 * p
    }
    n_bad <- 19  # number of outliers
    K <- 40  # number of blocks for MOM-CORE
    nloop <- 100  # number of replicates
    SNR <- 4  # signal-to-noise ratio
    num_method <- 6
    beta_true <- rep(1, p)  # true coefficients
    
    mse_meta <- array(0, dim = c(nloop, num_method, length(sub_meta)))
    mse_temp <- matrix(0, num_method, length(sub_meta))
    
    ##########################################################################
    for (i in 1:nloop) {
      
      ## Set the covariance matrix of X ##
      Sigma <- diag(p)
      for (ii in 1:p) {
        for (jj in 1:p) {
          Sigma[ii, jj] <- 0.6^(abs(ii - jj))
        }
      }
      set.seed(100 + 123 * i)
      
      ## Generate the design matrix X ##
      if (case[1] == "D1") {
        X <- mvrnorm(N, rep(0, p), Sigma)
      } else if (case[1] == "D2") {
        X <- mvrnorm(N, rep(0, p), Sigma)
        X <- exp(X)
        X <- X - matrix(rep(colMeans(X), N), N, p, byrow = TRUE)
      } else if (case[1] == "D3") {
        X <- mvtnorm::rmvt(N, Sigma, delta = rep(0, p), df = 3)
      }
      
      ## Sparsify the design matrix X ##
      if (case[2] == "R2") {
        X <- xSparse(X, 0.2)
      } else if (case[2] == "R3") {
        X <- xSparse(X, 0.4)
      } else if (case[2] == "R4") {
        X <- xSparse(X, 0.6)
      } else if (case[2] == "R5") {
        X <- xSparse(X, 0.8)
      }
      # X = X + matrix(runif(N * p, -0.1, 0.1), N, p)  # numerically sparse
      
      id_nz <- which(rowSums(abs(X)) != 0)
      X <- X[id_nz, ]
      Y <- yGenerate(X, beta_true, SNR)  # generate the response y
      N <- length(Y)
      
      ## Corrupt the generated data (or not) ##
      if (is_corrupt == TRUE) {
        data_corrupt <- dataCorrupt(n_bad, X, Y, beta_true, SNR)
        X <- data_corrupt[[1]]
        Y <- data_corrupt[[2]]
        
        ind_shuffle <- sample(N, N, replace = FALSE)
        X <- X[ind_shuffle, ]
        Y <- Y[ind_shuffle]
      }
      
      ## Compute the full sample OLS estimation ##
      svdf <- fast.svd(X)
      U <- svdf$u
      D <- svdf$d
      V <- svdf$v
      beta_ols <- as.vector(V %*% diag(1/D) %*% t(U) %*% Y)
      
      ## Compute the sampling probabilities for BLEV and SLEV ##
      PP <- rowSums(U^2)
      prob_blev <- PP/p
      prob_slev <- 0.9 * prob_blev + 0.1/N
      rm(svdf, PP, U, V, D)
      
      for (j in 1:length(sub_meta)) {
        set.seed(100 + 421 * j + 123 * i)
        sub <- sub_meta[j]
        
        ## Subsampling methods ##
        # UNIF
        id <- sample(N, sub, replace = FALSE)
        xt <- X[id, ]
        yt <- Y[id]
        fit <- lm(yt ~ xt - 1)
        unif_temp <- fit$coefficients
        rm(id, xt, yt, fit)
        
        # BLEV
        id <- sample(N, sub, replace = FALSE, prob_blev)
        xt <- X[id, ]
        yt <- Y[id]
        wgt <- 1/prob_blev[id]
        fit <- lm(yt ~ xt - 1, weights = wgt)
        blev_temp <- fit$coefficients
        rm(id, xt, yt, wgt, fit)
        
        # SLEV
        id <- sample(N, sub, replace = FALSE, prob_slev)
        xt <- X[id, ]
        yt <- Y[id]
        wgt <- 1/prob_slev[id]
        fit <- lm(yt ~ xt - 1, weights = wgt)
        slev_temp <- fit$coefficients
        rm(id, xt, yt, wgt, fit)
        
        # IBOSS
        id <- NULL
        X_c <- X
        ssub <- floor(sub/p/2)
        
        for (k in 1:p) {
          X_c[id, ] <- NA
          thres1 <- Rfast::nth(X_c[, k], ssub, descending = TRUE)
          id1 <- which(X_c[, k] >= thres1)[1:ssub]
          thres2 <- Rfast::nth(X_c[, k], ssub, descending = FALSE)
          id2 <- which(X_c[, k] <= thres2)[1:ssub]
          id <- c(id, id1, id2)
        }
        xt <- X[id, ]
        yt <- Y[id]
        
        fit <- lm(yt ~ xt - 1)
        iboss_temp <- fit$coefficients
        rm(X_c, id1, id2, id, xt, yt, fit)
        
        # CORE
        XX <- matrix(0, N, p)
        idd <- NULL
        
        thres <- colnth(abs(X), rep(sub, p), descending = TRUE, parallel = FALSE)
        for (kk in 1:p) {
          col <- abs(X[, kk])
          id <- which(col >= thres[kk])[1:sub]
          XX[id, kk] <- X[id, kk]
          idd <- unique(c(idd, id))
        }
        xxt <- XX[idd, ]
        yyt <- Y[idd]
        xt <- X[idd, ]
        
        xt1 <- Matrix(xt, sparse = TRUE)
        xxt1 <- Matrix(xxt, sparse = TRUE)
        core_temp <- eigenMultSolveSp(xt1, xxt1, yyt)
        rm(XX, col, thres, id, idd, xt, xxt, yyt, xt1, xxt1)
        
        # MOM-CORE
        core_dc_temp <- matrix(0, K, p)
        M <- floor(N/K)
        for (k in 1:K) {
          X_k <- X[((k - 1) * M + 1):(k * M), ]
          Y_k <- Y[((k - 1) * M + 1):(k * M)]
          XX_k <- matrix(0, M, p)
          idd <- NULL
          thres <- colnth(abs(X_k), rep(ceil(sub/K), p), descending = TRUE, parallel = FALSE)
          for (kk in 1:p) {
            col <- abs(X_k[, kk])
            id <- which(col >= thres[kk])[1:ceil(sub/K)]
            XX_k[id, kk] <- X_k[id, kk]
            idd <- unique(c(idd, id))
          }
          xxt_k <- XX_k[idd, ]
          yyt_k <- Y_k[idd]
          xt_k <- X_k[idd, ]
          
          xt1_k <- Matrix(xt_k, sparse = TRUE)
          xxt1_k <- Matrix(xxt_k, sparse = TRUE)
          core_dc_temp[k, ] <- eigenMultSolveSp(xt1_k, xxt1_k, yyt_k)
        }
        
        core_dc_temp <- na.omit(core_dc_temp)
        # core_dc = colMedians(core_dc_temp, na.rm = TRUE, parallel = FALSE)
        core_dc <- spat.med(core_dc_temp, tol = 1e-09)
        
        ######################################################
        coef <- cbind(unif_temp, blev_temp, slev_temp, iboss_temp,
                      core_temp, core_dc)
        coef[is.na(coef)] <- 0
        coef_mse <- colMeans((coef - beta_true)^2)/mean(beta_true^2)
        mse_temp[, j] <- coef_mse
      }
      
      mse_meta[i, , ] <- mse_temp
    }
    
    UNIF <- log(mse_meta[, 1, ])
    BLEV <- log(mse_meta[, 2, ])
    SLEV <- log(mse_meta[, 3, ])
    IBOSS <- log(mse_meta[, 4, ])
    CORE <- log(mse_meta[, 5, ])
    CORE_MOM <- log(mse_meta[, 6, ])
    
    ## Plot the results ##
    if (is_corrupt == TRUE) {
      mse_mat <- data.frame(mse = c(apply(UNIF, 2, mean), apply(BLEV, 2, mean), apply(SLEV, 2, mean), 
                                    apply(IBOSS, 2, mean), apply(CORE, 2, mean), apply(CORE_MOM, 2, mean)), 
                            sd = c(apply(UNIF, 2, sd), apply(BLEV, 2, sd), apply(SLEV, 2, sd), 
                                   apply(IBOSS, 2, sd), apply(CORE, 2, sd), apply(CORE_MOM, 2, sd)), 
                            Method = factor(rep(c("UNIF", "BLEV", "SLEV", "IBOSS", "CORE", "MOM-CORE"), 
                                                each = length(sub_meta))),
                            sub = rep(sub_meta/p, num_method))
      mse_mat$Method <- factor(mse_mat$Method, levels = c("UNIF", "BLEV", "SLEV", "IBOSS", "CORE", "MOM-CORE"))
      
      if (fig_num == 15) {
        leg_pos <- c(0.6, 0.69)
      } else {
        leg_pos <- ""
      }
      pd <- position_dodge(0.6)
      p1 <- ggplot(mse_mat, aes(x = sub, y = mse, group = Method, colour = Method))
      p23 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), 
                                     panel.border = element_rect(colour = "black"),
                                     axis.text = element_text(size = 19), 
                                     axis.title = element_text(size = 19),
                                     legend.position = leg_pos, 
                                     legend.title = element_blank(),
                                     legend.background = element_rect(fill = alpha("white", 0.6), color = "black"), 
                                     legend.key.width = unit(3, "line"),
                                     legend.key.height = unit(1, "line"), 
                                     legend.text = element_text(size = 13)) +
        geom_errorbar(aes(ymin = mse - sd, ymax = mse + sd), width = 1.5, position = pd, show.legend = FALSE) + 
        geom_line(aes(linetype = Method, size = Method), position = pd) + 
        geom_point(position = pd, aes(shape = Method), size = 2) + 
        scale_shape_manual(values = c(8, 4, 5, 7, 1, 2)) + 
        scale_linetype_manual(values = c(5, 4, 3, 2, 1, 1)) + 
        scale_size_manual(values = c(1, 1, 1, 1, 1, 1)) + 
        scale_color_manual(values = c("#999999", "#0033CC", "#3399FF", "orchid", "red", "#FF9900")) + 
        labs(title = paste(case, collapse = ", "), x = "r/p", y = "log(MSE)") + 
        theme(plot.title = element_text(hjust = 0.5, size = 19))
      
    } else {
      mse_mat <- data.frame(mse = c(apply(UNIF, 2, mean), apply(BLEV, 2, mean), apply(SLEV, 2, mean), 
                                    apply(IBOSS, 2, mean), apply(CORE, 2, mean)), 
                            sd = c(apply(UNIF, 2, sd), apply(BLEV, 2, sd), apply(SLEV, 2, sd), 
                                   apply(IBOSS, 2, sd), apply(CORE, 2, sd)), 
                            Method = factor(rep(c("UNIF", "BLEV", "SLEV", "IBOSS", "CORE"), 
                                                each = length(sub_meta))), 
                            sub = rep(sub_meta/p, num_method - 1))
      mse_mat$Method <- factor(mse_mat$Method, levels = c("UNIF", "BLEV", "SLEV", "IBOSS", "CORE"))
      
      if (fig_num == 15) {
        leg_pos <- c(0.73, 0.73)
      } else {
        leg_pos <- ""
      }
      pd <- position_dodge(0.4)
      p1 <- ggplot(mse_mat, aes(x = sub, y = mse, group = Method, colour = Method))
      p23 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), 
                                     panel.border = element_rect(colour = "black"),
                                     axis.text = element_text(size = 19), 
                                     axis.title = element_text(size = 19),
                                     legend.position = leg_pos, 
                                     legend.title = element_blank(),
                                     legend.background = element_rect(fill = alpha("white", 0.6), color = "black"), 
                                     legend.key.width = unit(3, "line"), 
                                     legend.key.height = unit(1, "line"), 
                                     legend.text = element_text(size = 13)) + 
        geom_errorbar(aes(ymin = mse - sd, ymax = mse + sd), width = 1, position = pd, show.legend = FALSE) +
        geom_line(aes(linetype = Method, size = Method), position = pd) +
        geom_point(position = pd, aes(shape = Method), size = 2) +
        scale_shape_manual(values = c(8, 4, 5, 7, 1)) + 
        scale_linetype_manual(values = c(5, 4, 3, 2, 1)) + 
        scale_size_manual(values = c(1, 1, 1, 1, 1)) + 
        scale_color_manual(values = c("#999999", "#0033CC", "#3399FF", "orchid", "red")) + 
        labs(title = paste(case, collapse = ", "), x = "r/p", y = "log(MSE)") + 
        theme(plot.title = element_text(hjust = 0.5, size = 19))
    }
    
    print(p23)
    fig_all[[fig_num]] <- p23
  }
}


pdf(file_name, width = 17.5, height = 10.5)
grid.arrange(fig_all[[1]], fig_all[[2]], fig_all[[3]], fig_all[[4]], fig_all[[5]],
             fig_all[[6]], fig_all[[7]], fig_all[[8]], fig_all[[9]], fig_all[[10]],
             fig_all[[11]], fig_all[[12]], fig_all[[13]], fig_all[[14]], fig_all[[15]],
             nrow = 3, ncol = 5)
dev.off()
