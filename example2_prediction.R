rm(list = ls())
gc()
## Import R Packages ##
list.of.packages <- c("data.table", "dplyr", "MASS", "corpcor", "pracma",
                      "Rfast", "Rcpp", "Matrix", "ggplot2", "gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)
library(data.table)
library(dplyr)
library(MASS)
library(corpcor)
library(pracma)
library(Rfast)
library(Rcpp)
library(Matrix)
library(ggplot2)
library(gridExtra)

## Import the source file ##
setwd(pwd())
sourceCpp("source/funcs.cpp")

## Load data ##
load("data/higgs.RData")
data <- data[, -1]
X_true <- as.matrix(data[, 1:21])  # observed design matrix
Y_true <- data[, 24]  # observed response
N <- nrow(X_true)  # sample size
p <- ncol(X_true)  # dimension

## Standardize data ##
X_bar <- colMeans(X_true)
X_sd <- apply(X_true, 2, sd)
X_std <- matrix(0, N, p)  # standardized design matrix
for (i in 1:N) {
  X_std[i, ] <- (X_true[i, ] - X_bar)/X_sd
}
Y_std <- (Y_true - mean(Y_true))/sd(Y_true)  # standardized response
rm(data, X_true, Y_true, X_bar, X_sd)

## Shuffle data ##
set.seed(2000)
ind_shuffle <- sample(N, N, replace = FALSE)
X_std <- X_std[ind_shuffle, ]
Y_std <- Y_std[ind_shuffle]

## Remove outliers ##
U <- fast.svd(X_std)$u
PP <- rowSums(U^2)
X_out <- boxplot(PP, plot = FALSE)$out
ind_out <- which(PP %in% X_out)
X_std <- X_std[-ind_out, ]
Y_std <- Y_std[-ind_out]
N <- nrow(X_std)
rm(ind_shuffle, U, PP, X_out, ind_out)

## Fit the model ##
fit <- lm(Y_std ~ X_std - 1)
beta_ols <- fit$coefficients
Y_pred <- X_std %*% beta_ols  # predicted response
Res <- Y_std - Y_pred  # residuals

## Set parameters ##
sub_meta <- c(2^4, 2^6, 2^8, 2^10, 2^12) * p  # subsample parameter r
nloop <- 100  # number of replicates
num_method <- 7
N_train <- floor(0.7 * N)  # training set size
N_test <- N - N_train  # testing set size
mse_meta <- array(0, dim = c(nloop, num_method, length(sub_meta)))
mse_temp <- matrix(0, num_method, length(sub_meta))

## Beginning ##
for(i in 1:nloop){
  set.seed(200 + 123 * i)
  ## Bootstrap ##
  id <- sample(1:N, N, replace = TRUE)  
  X <- X_std
  Y <- Y_pred + Res[id]
  
  ## Partition the training and testing sets ##
  id <- sample(1:N, N_train, replace = FALSE)
  X_test <- X[-id, ]
  Y_test <- Y[-id]
  X <- X[id, ]
  Y <- Y[id]
  
  ## Full sample OLS ##
  fit <- lm(Y ~ X - 1)
  ols_temp <- fit$coefficients
  rm(fit)

  ## Compute the sampling probabilities for BLEV and SLEV ##
  U <- fast.svd(X)$u
  PP <- rowSums(U^2)
  prob_blev <- PP/p
  prob_slev <- 0.9 * prob_blev + 0.1/nrow(X)
  rm(id, U, PP)
  
  for(j in 1:length(sub_meta)){
    print(c(i, j))
    set.seed(100 + 821 * j + 123 * i)
    sub <- sub_meta[j]
    
    ## Subsampling methods ## 
    # UNIF
    id <- sample(nrow(X), sub, replace = TRUE)
    xt <- X[id, ]
    yt <- Y[id]
    fit <- lm(yt ~ xt - 1)
    unif_temp <- fit$coefficients
    rm(id, xt, yt, fit)
    
    # BLEV
    id <- sample(nrow(X), sub, replace = TRUE, prob_blev)
    xt <- X[id, ]
    yt <- Y[id]
    wgt <- 1/prob_blev[id]
    fit <- lm(yt ~ xt - 1, weights = wgt)
    blev_temp <- fit$coefficients
    rm(id, xt, yt, fit, wgt)
    
    # SLEV
    id <- sample(nrow(X), sub, replace = TRUE, prob_slev)
    xt <- X[id, ]
    yt <- Y[id]
    wgt <- 1/prob_slev[id]
    fit <- lm(yt ~ xt - 1, weights = wgt)
    slev_temp <- fit$coefficients
    rm(id, xt, yt, fit, wgt)
    
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
    XX <- matrix(0, nrow(X), p)
    idd <- NULL
    thres <- colnth(abs(X), rep(sub, p), descending = TRUE, parallel = FALSE)
    for (kk in 1:p) {
      col <- abs(X[, kk])
      id1 <- which(col > thres[kk])
      id2 <- which(col == thres[kk])[1:(sub - length(id1))]
      id <- c(id1, id2)
      XX[id, kk] <- X[id, kk]
      idd <- unique(c(idd, id))
    }
    xxt <- XX[idd, ]
    yyt <- Y[idd]
    xt <- X[idd, ]
    
    xt1 <- Matrix(xt, sparse = TRUE)
    xxt1 <- Matrix(xxt, sparse = TRUE)
    core_temp <- eigenMultSolveSp(xt1, xxt1, yyt)
    rm(XX, col, thres, id1, id2, id, idd, xt, xxt, yyt, xt1, xxt1)
    
    # MOM-CORE
    K <- 5  # number of blocks
    core_dc_temp <- matrix(0, K, p)
    M <- floor(nrow(X)/K)
    for (k in 1:K) {
      X_k <- X[((k - 1) * M + 1):(k * M), ]
      Y_k <- Y[((k - 1) * M + 1):(k * M)]
      
      XX_k <- matrix(0, M, p)
      idd <- NULL
      thres <- colnth(abs(X_k), rep(ceil(sub/K), p), descending = TRUE, parallel = FALSE)
      for (kk in 1:p) {
        col <- abs(X_k[, kk])
        id1 <- which(col > thres[kk])
        id2 <- which(col == thres[kk])[1:(ceil(sub/K) - length(id1))]
        id <- c(id1, id2)
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
    core_dc <- spat.med(core_dc_temp, tol = 1e-09)
    rm(XX_k, X_k, Y_k, thres, id1, id2, id, idd, xt_k, xxt_k, yyt_k, xt1_k, xxt1_k)
    
    
    ######################################################
    coef <- cbind(ols_temp, unif_temp, blev_temp, slev_temp, iboss_temp, 
                  core_temp, core_dc)
    coef[is.na(coef)] <- 0
    coef_mse <- colMeans((X_test %*% coef - Y_test)^2)/mean(Y_test^2)
    mse_temp[, j] <- coef_mse
  }
  
  mse_meta[i, , ] <- mse_temp
}

## Plot the results ##
OLS <- log(mse_meta[, 1, ])
UNIF <- log(mse_meta[, 2, ])
BLEV <- log(mse_meta[, 3, ])
SLEV <- log(mse_meta[, 4, ])
IBOSS <- log(mse_meta[, 5, ])
CORE <- log(mse_meta[, 6, ])
CORE_MOM <- log(mse_meta[, 7, ])

mse_mat <- data.frame(mse = c(apply(OLS, 2, mean), apply(UNIF, 2, mean), apply(BLEV, 2, mean), 
                              apply(SLEV, 2, mean), apply(IBOSS, 2, mean), 
                              apply(CORE, 2, mean), apply(CORE_MOM, 2, mean)), 
                      sd = c(apply(OLS, 2, sd), apply(UNIF, 2, sd), apply(BLEV, 2, sd), 
                             apply(SLEV, 2, sd), apply(IBOSS, 2, sd), 
                             apply(CORE, 2, sd), apply(CORE_MOM, 2, sd)), 
                      Method = factor(rep(c("FullOLS", "UNIF", "BLEV", "SLEV", "IBOSS", "CORE", "MOM-CORE"), 
                                          each = length(sub_meta))),
                      sub = rep(log(sub_meta/p), num_method))
mse_mat$Method <- factor(mse_mat$Method, levels = c("FullOLS", "UNIF", "BLEV", "SLEV", "IBOSS", "CORE", "MOM-CORE"))

pd <- position_dodge(0.1)
p1 <- ggplot(mse_mat, aes(x = sub, y = mse, group = Method, colour = Method))
p23 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.border = element_rect(colour = "black"), 
                               axis.text = element_text(size = 17),
                               axis.title = element_text(size = 17), 
                               legend.position = "right", 
                               legend.title = element_blank(),
                               legend.key.width = unit(3, "line"), 
                               legend.key.height = unit(1.5, "line"),
                               legend.text = element_text(size = 15)) + 
  geom_errorbar(aes(ymin = mse - sd, ymax = mse + sd), width = 0.5, position = pd, show.legend = FALSE) +
  geom_line(aes(linetype = Method, size = Method), position = pd) + 
  geom_point(position = pd, aes(shape = Method), size = 2) + 
  scale_shape_manual(values = c(9, 8, 4, 5, 7, 1, 2)) + 
  scale_linetype_manual(values = c(6, 5, 4, 3, 2, 1, 1)) + 
  scale_size_manual(values = c(1, 1, 1, 1, 1, 1, 1)) + 
  scale_color_manual(values = c("grey2", "#999999", "#0033CC", "#3399FF", "orchid", "red", "#FF9900")) + 
  labs(x = "log(r/p)", y = "log(PMSE)") + 
  theme(plot.title = element_text(hjust = 0.5, size = 19))
p23
ggsave(filename = "example2_pmse.png", width = 6, height = 3.5, units = "in", dpi = 300)
