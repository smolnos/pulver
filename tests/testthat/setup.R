dir.create("data", showWarnings = FALSE)

test_that <- function(...) {
    testthat::test_that(...)
    file.remove(dir("data", full.names = TRUE))
}

expect_same_contents <- function(actual, expected, absolute_tolerance = 1e-8, order_matters = TRUE) {
    if (ncol(actual) != ncol(expected)) {
        fail(paste0("actual and expected have different numbers of columns: ",
            ncol(actual), " vs ", ncol(expected)))
        return()
    }
    if (!(if (order_matters) identical else setequal)(names(actual), names(expected))) {
        fail(paste0("different columns in actual and expected: (",
            paste0(names(expected), collapse = ", "), ") vs (", paste(names(actual), collapse = ", "), ")"))
        return()
    }
    if (nrow(actual) != nrow(expected)) {
        fail(paste0("actual and expected have different number of rows: ",
            nrow(actual), " vs ", nrow(expected)))
        return()
    }
    common_columns <- setdiff(names(expected), "pvalue")
    x <- merge(actual, expected, by = common_columns, suffixes = c(".actual", ".expected"))
    if (nrow(x) == 0) {
        fail(paste("Columns", paste(common_columns, collapse = ", "), "of actual and expected contain no common entries."))
        return()
    }
    if (nrow(x) != nrow(expected)) {
        fail(paste0("Duplicates in columns: ", paste(common_columns, collapse = ", ")))
        return()
    }
    for (col in common_columns) {
        if (any(missing <- !actual[[col]] %in% x[[col]])) {
            fail(paste0("Columns did not match ", paste0("(", paste(actual[[1]][missing], actual[[2]][missing], actual[[3]][missing], sep = ","), ")")))
            return()
        }
    }
    if (any(different <- abs(x$pvalue.actual - x$pvalue.expected) > absolute_tolerance)) {
        fail(paste0("actual and expected p-values differ:\n",
            paste(head(x$pvalue.actual[different]), head(x$pvalue.expected[different]), collapse = "\n", sep = " vs ")))
        return()
    }
    succeed("")
}

expect_setequal <- function(actual, expected) {
    if (!setequal(actual, expected)) {
        if (length(a_not_e <- setdiff(actual, expected)))
            fail(paste("In actual but not in expected:", paste(a_not_e, collapse = ", ")))
        if (length(e_not_a <- setdiff(expected, actual)))
            fail(paste("In expected but not in actual:", paste(e_not_a, collapse = ", ")))
        succeed()
    }
    succeed()
}

expect_no_error <- function(expr) {
    tryCatch(expr, error = function(c) {
        fail(paste("Got an error:", conditionMessage(c)))
    })
    succeed()
}
