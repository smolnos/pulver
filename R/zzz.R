#' @useDynLib pulver
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function(libpath) {
    library.dynam.unload("pulver", libpath)
}
