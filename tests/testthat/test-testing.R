source("setup.R")

context("create_xyz_files")

test_that("we get a list of 3 vectors with file names", {
    filenames <- create_xyz_files()
    expect_equal(length(filenames), 3)
    expect_is(filenames$x, "character")
})

test_that("we can get multiple files on demand", {
    filenames <- create_xyz_files(y = 1, x = 2, z = 3)
    expect_equal(sapply(filenames, length), c(y = 1, x = 2, z = 3))
})

test_that("there are no duplicated filenames", {
    filenames <- create_xyz_files()
    expect_true(!anyDuplicated(unlist(filenames)))
})

test_that("we produce txt files by default", {
    filenames <- create_xyz_files()
    expect_setequal(tools::file_ext(unlist(filenames)), "txt")
})

test_that("we can produce fvi files on demand", {
    filenames <- create_xyz_files(type = "fvi")
    expect_setequal(tools::file_ext(unlist(filenames)), "fvi")
})

test_that("we can produce different types of files at the same time", {
    filenames <- unlist(create_xyz_files(type = c("txt", "fvi")))
    expect_true("txt" %in% tools::file_ext(filenames))
    expect_true("fvi" %in% tools::file_ext(filenames))
})

test_that("we can request a specific number of columns", {
    filenames <- create_xyz_files(ycols = 1, xcols = 2, zcols = 3)
    expect_equal(ncol(read_data(filenames$y)), 1)
    expect_equal(ncol(read_data(filenames$x)), 2)
    expect_equal(ncol(read_data(filenames$z)), 3)
})

test_that("a given set of files contains no duplicate column names", {
    matrices <- lapply(create_xyz_files(y = 2)$y, read.table)
    column_names <- unlist(lapply(matrices, colnames))
    expect_true(!anyDuplicated(column_names))
})

test_that("x, y, and z files use different column names", {
    filenames <- create_xyz_files()
    ynames <- names(read.table(filenames$y))
    xnames <- names(read.table(filenames$x))
    expect_true(length(intersect(xnames, ynames)) == 0)
})

test_that("by default all files are saved in data/ directory", {
    filenames <- create_xyz_files()
    expect_setequal(dirname(unlist(filenames)), "data")
})

test_that("we can set the number of rows", {
    filenames <- create_xyz_files(nobs = 5)
    expect_equal(nrow(read_data(filenames$y)), 5)
})

context("write_datasets")

test_that("we can supply multiple files", {
    write_datasets(c("data/y", "data/x"))
    expect_true(file.exists("data/y"))
    expect_true(file.exists("data/x"))
})

test_that("we can create fvi files", {
    write_datasets("data/y.fvi", ncol = 2)
    expect_equal(ncol(read_data("data/y.fvi")), 2)
})

test_that("we can set the number of rows", {
    write_datasets("data/y.txt", nobs = 5)
    expect_equal(nrow(read.delim("data/y.txt")), 5)
})

context("create_xyz_matrices")

test_that("happy path", {
    data <- create_xyz_matrices(ycols = 2, xcols = 3, zcols = 4, nobs = 10)
    expect_equal(ncol(data$y), 2)
    expect_equal(ncol(data$x), 3)
    expect_equal(ncol(data$z), 4)
    expect_equal(nrow(data$y), 10)
    expect_equal(nrow(data$x), 10)
    expect_equal(nrow(data$z), 10)
})

context("write_databel and read_databel")

test_that("we can write and read DatABEL files", {
    data <- matrix(as.double(1:12), 3)
    write_databel(data, "data/file")
    expect_equal(read_databel("data/file.fvi"), data)
})
