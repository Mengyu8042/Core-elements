yGenerate <- function(x, beta, SNR, noise_type = NULL) {
  # generate the response vector
  y_raw <- as.vector(x %*% beta)
  N <- length(y_raw)
  if (is.null(noise_type)) {
    h_term <- 0
  } else {
    if (noise_type == "M1") {
      h_term <- x[, 3] * x[, 8]
    } else if (noise_type == "M2") {
      h_term <- x[, 3] * sin(x[, 8])
    } else if (noise_type == "M3") {
      h_term <- x[, 3]^2
    } 
    h_term <- h_term/max(abs(h_term))
  }
  y <- y_raw + rnorm(N, 10 * h_term, sd(y_raw)/sqrt(SNR))
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

dataCorrupt <- function(n_bad, x, y, beta) {
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
  
  return(list(X = x_corrupt, Y = y_corrupt))
}

dataGenerate <- function(N, p, Sigma, beta_true, SNR, case = c("D1", "R1"), 
                         is_corrupt = FALSE, n_bad = 0,
                         noise_type = NULL) {
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
  X <- X + matrix(runif(N * p, -0.1, 0.1), N, p)  # numerically sparse
  
  id_nz <- which(rowSums(abs(X)) != 0)
  X <- X[id_nz, ]
  Y <- yGenerate(X, beta_true, SNR, noise_type)  # generate the response y
  
  ## Corrupt the generated data (or not) ##
  if (is_corrupt) {
    data_corrupt <- dataCorrupt(n_bad, X, Y, beta_true)
    X <- data_corrupt$X
    Y <- data_corrupt$Y
    
    ind_shuffle <- sample(N, N, replace = FALSE)
    X <- X[ind_shuffle, ]
    Y <- Y[ind_shuffle]
  }
  
  return(list(X = X, Y = Y))
}

subsample_lm <- function(X, Y, id, weights = NULL) {
  xt <- X[id, ]
  yt <- Y[id]
  fit <- lm(yt ~ xt - 1, weights = weights)
  return(fit$coefficients)
}

method_UNIF <- function(X, Y, sub) {
  id <- sample(nrow(X), sub, replace = FALSE)
  beta_est <- subsample_lm(X, Y, id)
  return(beta_est)
}

method_BLEV <- function(X, Y, sub, prob_blev) {
  id <- sample(nrow(X), sub, replace = FALSE, prob_blev)
  wgt <- 1/prob_blev[id]
  beta_est <- subsample_lm(X, Y, id, weights = wgt)
  return(beta_est)
}

method_SLEV <- function(X, Y, sub, prob_slev) {
  id <- sample(nrow(X), sub, replace = FALSE, prob_slev)
  wgt <- 1/prob_slev[id]
  beta_est <- subsample_lm(X, Y, id, weights = wgt)
  return(beta_est)
}

method_IBOSS <- function(X, Y, sub) {
  N <- nrow(X)
  p <- ncol(X)
  id <- NULL
  X_c <- X
  ssub <- ceil(sub/p/2)
  
  for (k in 1:p) {
    X_c[id, ] <- NA
    thres1 <- Rfast::nth(X_c[, k], ssub, descending = TRUE, na.rm = TRUE)
    id1 <- which(X_c[, k] >= thres1)[1:ssub]
    thres2 <- Rfast::nth(X_c[, k], ssub, descending = FALSE, na.rm = TRUE)
    id2 <- which(X_c[, k] <= thres2)[1:ssub]
    id <- c(id, id1, id2)
  }
  if (length(id) > sub) {
    id <- sample(id, sub, replace = FALSE)
  }
  beta_est <- subsample_lm(X, Y, id)
  return(beta_est)
}

method_OSS <- function(X, Y, sub) {
  id <- OAJ2_cpp(X, sub)
  beta_est <- subsample_lm(X, Y, id)
  return(beta_est)
}

method_DOUBLY <- function(X, Y, sub) {
  sub2 <- 5 * sub
  id <- sample(nrow(X), sub2, replace = TRUE)
  xt <- X[id, ]
  yt <- Y[id]
  S_cw <- zeros(sub, sub2)
  for (jj in 1:ncol(S_cw)) {
    S_cw[sample(sub, 1), jj] <- sample(c(1, -1), 1)
  }
  S_cw_sp <- Matrix(S_cw, sparse = TRUE)
  xt <- as.matrix(S_cw_sp %*% xt)
  yt <- as.matrix(S_cw_sp %*% yt)
  beta_est <- lm(yt ~ xt - 1)$coefficients
  return(beta_est)
}

method_DOPT <- function(X, Y, sub, ord_dopt) {
  id <- ord_dopt[1:sub]
  beta_est <- subsample_lm(X, Y, id)
  return(beta_est)
}

method_CORE <- function(X, Y, sub) {
  N <- nrow(X)
  p <- ncol(X)
  XX <- matrix(0, N, p)
  thres <- colnth(abs(X), rep(sub, p), descending = TRUE, parallel = FALSE)
  idd <- NULL
  
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
  beta_est <- eigenMultSolveSp(xt1, xxt1, yyt)
  return(beta_est)
}

method_FULL <- function(X, Y) {
  svdf <- fast.svd(X)
  U <- svdf$u
  D <- svdf$d
  V <- svdf$v
  beta_est <- as.vector(V %*% diag(1/D) %*% t(U) %*% Y)
  
  PP <- rowSums(U^2)
  prob_blev <- PP/ncol(X)
  prob_slev <- 0.9 * prob_blev + 0.1/nrow(X)
  
  Z <- X - colMeans(X)
  U <- fast.svd(Z)$u
  PP <- rowSums(U^2)
  ord_dopt <- Order(PP, descending = TRUE)
  
  return(list(beta_est = beta_est, prob_blev = prob_blev,
              prob_slev = prob_slev, ord_dopt = ord_dopt))
}

plot_mse <- function(mse_mat, sub_meta, is_predict, method_colors, fig_title) {
  y_axis <- ifelse(is_predict, expression(log[10](PMSE)), expression(log[10](MSE)))
  
  width <- (max(sub_meta) - min(sub_meta))/p
  pd <- position_dodge(width/30)
  
  p1 <- ggplot(mse_mat, aes(x = sub, y = mse, group = Method, colour = Method))
  p23 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), 
                                 panel.border = element_rect(colour = "black"),
                                 axis.text = element_text(size = 19), 
                                 axis.title = element_text(size = 19),
                                 legend.position = "right", 
                                 legend.justification = c(1, 1), 
                                 legend.box.just = "right",
                                 legend.margin = margin(1, 1, 1, 1),
                                 legend.box.margin = margin(1, 1, 1, 1),
                                 legend.title = element_blank(),
                                 legend.background = element_rect(fill = alpha("white", 0.6), color = "black"), 
                                 legend.key.width = unit(3, "line"), 
                                 legend.key.height = unit(1, "line"), 
                                 legend.text = element_text(size = 17)) + 
    geom_errorbar(aes(ymin = mse - sd, ymax = mse + sd), width = 1, position = pd, show.legend = FALSE) +
    geom_line(aes(linetype = Method, size = Method), position = pd) +
    geom_point(position = pd, aes(shape = Method), size = 2) +
    scale_shape_manual(values = rev(1:length(method_colors))) + 
    scale_linetype_manual(values = rev(1:length(method_colors))) + 
    scale_size_manual(values = rep(1, length(method_colors))) + 
    scale_color_manual(values = method_colors) + 
    labs(title = fig_title, x = "r/p", y = y_axis) + 
    theme(plot.title = element_text(hjust = 0.5, size = 19))
  
  return(p23)
}

get_legend <- function(myplot) {
  tmp <- ggplotGrob(myplot + theme(legend.position = "bottom",
                                   legend.justification = "center",
                                   legend.box.just = "center") +
                      guides(color = guide_legend(nrow = 1, byrow = TRUE)))
  legend <- gtable::gtable_filter(tmp, "guide-box")
  return(legend)
}