Sys.setenv("R_TESTS" = "")
library(testthat)
library(pulver)

test_check("pulver")
