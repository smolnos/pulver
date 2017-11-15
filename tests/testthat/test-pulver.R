source("setup.R")
source("testing.R")

context("pulverize")

test_that("happy path", {
    mat <- create_xyz_matrices()
    actual <- pulverize(mat$y, mat$x, mat$z)
    expected <- pulverize_lm(mat$y, mat$x, mat$z)
    expect_same_contents(actual, expected)
})

test_that("we can suppress the return value", {
    mat <- create_xyz_matrices()
    expect_null(pulverize(mat$y, mat$x, mat$z, "data/out.txt", suppress_return = TRUE))
    actual <- utils::read.delim("data/out.txt", as.is = TRUE)
    expected <- pulverize_lm(mat$y, mat$x, mat$z)
    expect_same_contents(actual, expected)
})

test_that("we can write to an output file", {
    mat <- create_xyz_matrices()
    expected <- pulverize(mat$y, mat$x, mat$z, "data/out.txt")
    actual <- utils::read.delim("data/out.txt", as.is = TRUE)
    expect_same_contents(actual, expected)
})

test_that("we can filter by p-value", {
    mat <- create_xyz_matrices()
    actual <- pulverize(mat$y, mat$x, mat$z, pvalue_threshold = .5)
    expected <- pulverize_lm(mat$y, mat$x, mat$z, pvalue_threshold = .5)
    expect_same_contents(actual, expected)
})

test_that("we get an error if an output file exists already", {
    mat <- create_xyz_matrices()
    file.create("data/exists_already.txt")
    expect_error(pulverize(mat$y, mat$x, mat$z, "data/exists_already.txt"),
        "File exists already: data/exists_already.txt")
})

test_that("we get an error if one matrix does not have column names", {
    mat <- create_xyz_matrices()
    mat$x <- matrix(rnorm(10 * 3), ncol = 3)
    expect_error(pulverize(mat$y, mat$x, mat$z),
        "Column names for at least one matrix is missing.")
})

test_that("we can overwrite the output file on request", {
    mat <- create_xyz_matrices()
    file.create("data/out.txt")
    expect_no_error(pulverize(mat$y, mat$x, mat$z, "data/out.txt", overwrite = TRUE))
})

test_that("we can choose column names for the results table", {
    mat <- create_xyz_matrices()
    actual <- pulverize(mat$y, mat$x, mat$z, colnames = c("a", "b", "c"))
    expect_equal(names(actual), c("a", "b", "c", "pvalue"))
    expected <- pulverize_lm(mat$y, mat$x, mat$z)
    colnames(expected) <- c("a", "b", "c", "pvalue")
    expect_same_contents(actual, expected)
})

test_that("we can use multiple cores in parallel", {
    mat <- create_xyz_matrices()
    actual <- pulverize(mat$y, mat$x, mat$z, cores = 4)
    expected <- pulverize_lm(mat$y, mat$x, mat$z)
    expect_same_contents(actual, expected)
})

context("pulverize_all")

test_that("we can write to an output file", {
    mat <- create_multiple_xyz_matrices()
    expected <- pulverize_all(mat$y, mat$x, mat$z, "data/out.txt", suppress_return = FALSE, overwrite = TRUE)
    actual <- utils::read.delim("data/out.txt", as.is = TRUE)
    expect_same_contents(actual, expected)
})

test_that("we can select column names for the result table", {
    mat <- create_multiple_xyz_matrices()
    actual <- pulverize_all(mat$y, mat$x, mat$z, colnames = c("a", "b", "c"), suppress_return = FALSE)
    expected <- pulverize_lm_all(mat$y, mat$x, mat$z)
    names(expected) <- c("a", "b", "c", "pvalue")
    expect_same_contents(actual, expected)
})

test_that("we can supply multiple x, y, and z matrices", {
    mat <- create_multiple_xyz_matrices(y = 2, x = 2, z = 2)
    actual <- pulverize_all(mat$y, mat$x, mat$z, suppress_return = FALSE)
    expected <- pulverize_lm_all(mat$y, mat$x, mat$z)
    expect_same_contents(actual, expected)
})

test_that("we can supply a mixture of txt and fvi files", {
    files <- create_xyz_files(y = 2, type = c("txt", "fvi"))
    expect_true("txt" %in% tools::file_ext(files$y))
    expect_true("fvi" %in% tools::file_ext(files$y))
    actual <- pulverize_all(files$y, files$x, files$z, suppress_return = FALSE)
    expected <- pulverize_lm_all(files$y, files$x, files$z)
    expect_same_contents(actual, expected)
})

test_that("we get an error if an output file exists already", {
    mat <- create_multiple_xyz_matrices()
    file.create("data/exists_already.txt")
    expect_error(pulverize_all(mat$y, mat$x, mat$z, "data/exists_already.txt"),
        "File exists already: data/exists_already.txt")
})

test_that("we can overwrite the output file on request", {
    mat <- create_multiple_xyz_matrices()
    file.create("data/out.txt")
    expect_no_error(pulverize_all(mat$y, mat$x, mat$z, "data/out.txt", overwrite = TRUE))
})

context("pulverize_lm_all")

test_that("we can supply multiple x, y, z matrices", {
    mat <- create_multiple_xyz_matrices(y = 2, x = 2, z = 2)
    actual <- pulverize_lm_all(mat$y, mat$x, mat$z)
    index <- expand.grid(y = 1:2, x = 1:2, z = 1:2)
    expected <- do.call(rbind, Map(function(y, x, z) {
        pulverize_lm(mat$y[[y]], mat$x[[x]], mat$z[[z]])
    }, index$y, index$x, index$z))
    expect_same_contents(actual, expected)
})

test_that("we can select the column names of the result table", {
    mat <- create_multiple_xyz_matrices()
    actual <- pulverize_lm_all(mat$y, mat$x, rep(mat$z), colnames = c("a", "b", "c"))
    expect_equal(colnames(actual), c("a", "b", "c", "pvalue"))
})

context("r2t and t2r")

test_that("r2t is the inverse of t2r and vice versa", {
    expect_equal(t2r(r2t(+.3, 10), 10), +.3)
    expect_equal(r2t(t2r(-.3, 10), 10), -.3)
})

context("r2p and p2r")

test_that("r2p is the inverse of p2r and vice versa", {
    expect_equal(p2r(r2p(-.3, 10), 10), +.3)
    expect_equal(r2p(p2r(+.3, 10), 10), +.3)
})

context("create_pulverize_annotation")

test_that("create_pulverize_annotation", {
    y <- matrix(ncol = 2, dimnames = list(1, c("y1", "y2")))
    x <- matrix(ncol = 3, dimnames = list(1, c("x1", "x2", "x3")))
    z <- matrix(ncol = 2, dimnames = list(1, c("z1", "z2")))
    actual <- create_pulverize_annotation(y, x, z)
    expected <- as.data.frame(rbind(
        c("y1", "x1", "z1"),
        c("y1", "x1", "z2"),
        c("y1", "x2", "z1"),
        c("y1", "x2", "z2"),
        c("y1", "x3", "z1"),
        c("y1", "x3", "z2"),
        c("y2", "x1", "z1"),
        c("y2", "x1", "z2"),
        c("y2", "x2", "z1"),
        c("y2", "x2", "z2"),
        c("y2", "x3", "z1"),
        c("y2", "x3", "z2")),
        stringsAsFactors = FALSE)
    names(expected) <- c("y", "x", "z")
    expect_equal(actual, expected)
})

context("replace_NA_by_mean")
test_that("replace_NA_by_mean", {
    data <- create_xyz_matrices(ycols = 2, xcols = 3, zcols = 4, nobs = 10)
    data$x[10, 2] <- NA
    data$x[5, 3] <- NA
    actual <- replace_NA_by_mean(data$x)
    expected <- data$x
    expected[10, 2] <- mean(expected[, 2], na.rm = TRUE)
    expected[5, 3] <- mean(expected[, 3], na.rm = TRUE)
    expect_equal(actual, expected)
})

context("matrices_have_same_rownames")
test_that("matrices_have_rownames", {
    mat <- create_xyz_matrices(ycols = 2, xcols = 3, zcols = 4, nobs = 10)
    rownames(mat$x) <- NULL
    expect_warning(pulverize(mat$y, mat$x, mat$z), "Matrix xmat does not have row names.")
    rownames(mat$z) <- NULL
    expect_warning(pulverize(mat$y, mat$x, mat$z), "Matrices xmat, zmat do not have row names.")
})

test_that("matrices_have_same_rownames", {
    mat <- create_xyz_matrices(ycols = 2, xcols = 3, zcols = 4, nobs = 10)
    rownames(mat$x) <- 2:11
    expect_error(pulverize(mat$y, mat$x, mat$z),
        "Matrices do not have identical row names.")
})
