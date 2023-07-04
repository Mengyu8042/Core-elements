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

## Beginning ##
for (case_type in c("D1", "D2", "D3")) {
  for (spar_type in c("R1", "R2", "R3", "R4", "R5")) {
    fig_num <- fig_num + 1
    case <- c(case_type, spar_type)
    print(c(fig_num, case))
    
    ## Set parameters ##
    N <- 10000  # sample size n
    p <- 100
    sub_meta <- N - 100 * (1:10)
    nloop <- 100  # number of replicates
    SNR <- 4  # signal-to-noise ratio
    num_method <- 2
    beta_true <- rep(1, p)  # true coefficients

    epsY_meta <- array(0, dim = c(nloop, num_method, length(sub_meta)))
    epsX_meta <- matrix(0, nloop, length(sub_meta))

    ##########################################################################
    for (i in 1:nloop) {
      set.seed(100 + 123 * i)
      
      Sigma <- diag(p)
      if (case[1] == "D1") {
        X <- mvrnorm(N, rep(0, p), Sigma)
      } else if (case[1] == "D2") {
        X <- mvrnorm(N, rep(0, p), Sigma)
        X <- exp(X)
        X <- X - matrix(rep(colMeans(X), N), N, p, byrow = TRUE)
      } else if (case[1] == "D3") {
        X <- mvtnorm::rmvt(N, Sigma, delta = rep(0, p), df = 3)
      }
      
      if (case[2] == "R2") {
        X <- xSparse(X, 0.2)
      } else if (case[2] == "R3") {
        X <- xSparse(X, 0.4)
      } else if (case[2] == "R4") {
        X <- xSparse(X, 0.6)
      } else if (case[2] == "R5") {
        X <- xSparse(X, 0.8)
      }
      X <- X + matrix(runif(N * p, -0.1, 0.1), N, p)  # numerically sparse
  
      X_scale <- matrix(rep(apply(X, 2, sd), N), N, p, byrow = TRUE)
      X <- X/X_scale
      
      id_nz <- which(rowSums(abs(X)) != 0)
      X <- X[id_nz, ]
      Y <- yGenerate(X, beta_true, SNR)
      N <- length(Y)
      
      svdf <- fast.svd(X)
      D <- svdf$d
      lam_max <- D[1]
      lam_min <- D[p]
      kappa <- lam_max/lam_min
      
      fit <- lm(Y ~ X - 1)
      beta_ols <- fit$coefficients
      sse <- sum((Y - X %*% beta_ols)^2)
      Y_normSq <- sum(Y^2)
      err_normSq <- sum((Y - X %*% beta_true)^2)
      
      for(j in 1:length(sub_meta)){
        set.seed(100 + 421 * j + 123 * i)
        sub <- sub_meta[j]
        
        # Core-elements
        XX <- matrix(0, N, p)
        idd <- NULL
        
        thres <- colnth(abs(X), rep(sub, p), descending = TRUE, parallel = FALSE)
        for(kk in 1:p){
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
        rm(col, thres, id, idd, xt, xxt, yyt, xt1, xxt1)
        
        sse_core <- sum((Y - X %*% core_temp)^2)
        
        epsX_meta[i, j] <- norm(X - XX, type = "2")/lam_max
        
        epsY_prac <- sse_core/sse - 1
        
        frac1 <- (kappa^2 + 1)^2 * err_normSq
        frac2 <- (1/(epsX_meta[i, j] * kappa^2) - 1)^2 * sse
        epsY_theo <- frac1/frac2
        
        epsY_meta[i, , j] <- cbind(epsY_prac, epsY_theo)
      }
    }

    PRAC <- apply(log(epsY_meta[, 1, ]), 2, mean)
    THEO <- apply(log(epsY_meta[, 2, ]), 2, mean)
    
    ## Plot the results ##
    eps_mat <- data.frame(epsY = c(PRAC, THEO),
                          Method = factor(rep(c("Empirical", "Theoretical"), 
                                              each = length(sub_meta))),
                          epsX = rep(apply(log(epsX_meta), 2, mean), num_method))
    
    if (fig_num == 15) {
      leg_pos <- c(0.59, 0.85)
    } else {
      leg_pos <- ""
    }
    pd <- position_dodge(0)
    p1 <- ggplot(eps_mat, aes(x = epsX, y = epsY, group = Method, colour = Method))
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
                                   legend.text = element_text(size = 17)) +
      geom_line(aes(linetype = Method, size = Method), position = pd) +
      geom_point(position = pd, aes(shape = Method), size = 2) +
      scale_linetype_manual(values = c(1, 2)) +
      scale_size_manual(values = c(1, 1)) +
      scale_color_manual(values = c("red", "#0033CC")) +
      labs(title = paste(case, collapse = ", "), x = expression("log(ε')"), 
           y = expression("log(ε)")) +
      theme(plot.title = element_text(hjust = 0.5, size = 19))
    print(p23)
    fig_all[[fig_num]] <- p23
  }
}

pdf("error_bound.pdf", width = 17.5, height = 10.5)
grid.arrange(fig_all[[1]], fig_all[[2]], fig_all[[3]], fig_all[[4]], fig_all[[5]],
             fig_all[[6]], fig_all[[7]], fig_all[[8]], fig_all[[9]], fig_all[[10]],
             fig_all[[11]], fig_all[[12]], fig_all[[13]], fig_all[[14]], fig_all[[15]],
             nrow = 3, ncol = 5)
dev.off()