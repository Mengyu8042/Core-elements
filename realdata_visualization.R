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

## Load data ##
setwd(pwd())
load("data/scRNA_top3000.RData")
data <- as.matrix(data)
id_y <- order(colSums(data != 0), decreasing = TRUE)[1]
Y_true <- data[, id_y]  # observed response
data <- data[, -id_y]
X_true <- data[, order(apply(data, 2, sd), decreasing = TRUE)[1:500]]  # observed design matrix
id <- which(rowSums(X_true) != 0)
X_true <- X_true[id, ]
Y_true <- Y_true[id]
N <- nrow(X_true)  # sample size
p <- ncol(X_true)  # dimension

## Standardize data ##
X_scale <- matrix(rep(apply(X_true, 2, sd), N), N, p, byrow = TRUE)
X_std <- X_true/X_scale  # standardized design matrix
Y_std <- Y_true/sd(Y_true)  # standardized response
rm(data, X_true, Y_true, X_scale)

## Shuffle data ##
ind_shuffle <- sample(N, N, replace = FALSE)
X_std <- X_std[ind_shuffle, ]
Y_std <- Y_std[ind_shuffle]
rm(ind_shuffle)

## Visualize data ##
# Box plot - x
U <- fast.svd(X_std)$u
PP <- log(rowSums(U^2))
data_x <- data.frame(x = PP)
p1 <- ggplot(data_x, aes(y = x)) + geom_boxplot()
p2 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              axis.text=element_text(size = 18),
                              axis.title=element_text(size = 18)) +
  labs(x = "x", y = "log(leverage value)") +
  scale_x_discrete()

# Box plot - y
data_y <- data.frame(x = log(Y_std + 1))
p3 <- ggplot(data_y, aes(y = x)) + geom_boxplot()
p4 <- p3 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              axis.text = element_text(size = 18),
                              axis.title = element_text(size = 18)) +
  labs(x = "y", y = "log(response value)") +
  scale_x_discrete()

pdf("scRNA_boxplot.pdf", width = 5, height = 3.5)
grid.arrange(p2, p4, nrow = 1, ncol = 2)
dev.off()

# Histogram
X_std[X_std > quantile(X_std, 0.999)] <- NA
data_x <- data.frame(x = as.vector(X_std))
pdf("scRNA_hist.pdf", width = 3.5, height = 3.5)
p6 <- ggplot(data_x, aes(x = x))
p7 <- p6 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              axis.text = element_text(size = 17),
                              axis.title = element_text(size = 17),
                              legend.position = "none") +
  geom_histogram(aes(y = ..count../(N * p - sum(is.na(X_std))))) +
  labs(x = "value", y = "frequency") +
  xlim(c(-0.5, 11))
p7
dev.off()