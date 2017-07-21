pulverize_mixabel <- function(xmat, ymat, zmat, pvalue_threshold = 1, colnames = c("x", "y", "z")) {
    result <- cbind(pulver:::create_pulverize_annotation(xmat, zmat, ymat, colnames = colnames[c(1, 3, 2)]), t = NA_real_)
    threshold <- pulver:::p2t(pvalue_threshold, df = nrow(xmat) - 4)
    data <- as.data.frame(cbind(xmat, zmat))
    for (i in colnames(xmat)) {
        for (k in colnames(zmat)) {
            formula <- as.formula(with(data, sprintf("%s ~ SNP * %s", i, k)))
            Zx <- as.data.frame(MixABEL::GWFGLS(formula, data, genodata = ymat, include.means = FALSE))
            result[which(result$x == i & result$z == k), "t"] <- Zx[, 6]/Zx[, 8]
        }
    }
    result <- result[which(abs(result$t) > threshold),]
}

library(testthat)
library(pulver)
library(MixABEL, quietly = TRUE)

source("../setup.R")

context("pulverize_mixabel")

test_that("we can filter by p-value", {
    mat <- pulver:::create_xyz_matrices()
    actual <- pulverize_mixabel(mat$x, mat$y, mat$z, pvalue_threshold = .5)
    actual$pvalue <- pulver:::t2p(actual$t, df = nrow(mat$x) - 4)
    actual$t <- NULL
    expected <- pulver:::pulverize_lm(mat$x, mat$z, mat$y, pvalue_threshold = .5)
    colnames(expected) = c("x", "z", "y", "pvalue")
    expect_same_contents(actual, expected)
})
