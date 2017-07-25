## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----Github, results = "hide"--------------------------------------------
library(devtools)
install_github("smolnos/pulver")

## ----create_matrices-----------------------------------------------------
set.seed(369)
nobs <- 100
Y <- matrix(rnorm(nobs * 2), ncol = 2, dimnames = list(paste0("row", 1:nobs),
paste0("column", 1:2)))
X <- matrix(rnorm(nobs * 3), ncol = 3, dimnames = list(paste0("row", 1:nobs),
paste0("column", 1:3)))
Z <- matrix(rnorm(nobs * 4), ncol = 4, dimnames = list(paste0("row", 1:nobs),
paste0("column", 1:4)))

## ----pulver--------------------------------------------------------------
library(pulver)
pulverize(Y, X, Z)

## ----bonf----------------------------------------------------------------
pulverize(Y, X, Z, pvalue_threshold = 0.05/24)

## ----cores---------------------------------------------------------------
nobs <- 1000
Y <- matrix(rnorm(nobs * 20), ncol = 20, dimnames = list(paste0("row", 1:nobs),
paste0("column", 1:20)))
X <- matrix(rnorm(nobs * 300), ncol = 300, dimnames = list(paste0("row", 1:nobs),
paste0("column", 1:300)))
Z <- matrix(rnorm(nobs * 40), ncol = 40, dimnames = list(paste0("row", 1:nobs),
paste0("column", 1:40)))
system.time(pulverize(Y, X, Z))
system.time(pulverize(Y, X, Z, cores = 2))

## ----create_lists--------------------------------------------------------
nobs <- 100
Y <- matrix(rnorm(nobs * 20), ncol = 20, dimnames = list(NULL,
paste0("column", 1:20)))
Y1 <- Y[,1:10]
Y2 <- Y[,11:20]
write.table(Y1, file = "Y1.txt", sep = "\t")
write.table(Y2, file = "Y2.txt", sep = "\t")
Ylist <- as.list(c("Y1.txt", "Y2.txt"))

## ----databel, results = "hide"-------------------------------------------
X <- matrix(rnorm(nobs * 30), ncol = 30, dimnames = list(paste0("row", 1:nobs),
paste0("column", 1:30)))
X1 <- X[,1:15]
X2 <- X[,16:30]
write_databel(X1, "X1")
write_databel(X2, "X2")

## ----Xlist---------------------------------------------------------------
Xlist <- as.list(c("X1.fvi", "X2.fvi"))

## ----Zlist---------------------------------------------------------------
Z <- matrix(rnorm(nobs * 4), ncol = 4, dimnames = list(paste0("row", 1:nobs),
paste0("column", 1:4)))
Zlist <- list(NULL)
Zlist[[1]] <- Z

## ----pulverize_all-------------------------------------------------------
pulverize_all(Ylist, Xlist, Zlist, output_file = "output_file.txt")
head(read.delim("output_file.txt"))

## ----remove, results = "hide"--------------------------------------------
file.remove("Y1.txt", "Y2.txt", "X1.fvi" ,"X1.fvd", "X2.fvi", "X2.fvd", "output_file.txt")

