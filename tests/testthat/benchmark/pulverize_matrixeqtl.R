pulverize_matrixeqtl <- function(xmat, ymat, zmat, output_file, pvalue_threshold = 1) {
    # Create 3 SlicedData objects for the analysis
    y = SlicedData$new( t(ymat) )
    x = SlicedData$new( t(xmat) )
    zmat = t(zmat)
    for(zncol in seq_len(nrow(zmat))) {
        z = SlicedData$new(zmat[zncol, , drop = FALSE])
        # Call the main analysis function
        # result[[zcolname]]
        Matrix_eQTL_main(
            snps = y,
            gene = x,
            cvrt = z,
            output_file,
            pvOutputThreshold = pvalue_threshold,
            useModel = modelLINEAR_CROSS,
            verbose = FALSE,
            pvalue.hist = FALSE,
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = TRUE)
    }
}


library(testthat)
library(pulver)
library(MatrixEQTL)

source("../setup.R")

context("pulverize_eqtl")

test_that("we can filter by p-value", {
    mat <- pulver:::create_xyz_matrices(zcols = 1)
    output_file <- "result"
    pulverize_matrixeqtl(mat$x, mat$y, mat$z, output_file, pvalue_threshold = 0.75)
    actual <- read.delim("result")
    actual <- actual[, c("SNP", "gene", "p.value")]
    actual$z <- colnames(mat$z)
    names(actual) <- c("y", "x", "pvalue", "z")
    expected <- pulver:::pulverize_lm(mat$x, mat$y, mat$z, pvalue_threshold = 0.75)
    expect_same_contents(actual, expected, order_matters = FALSE)
})
