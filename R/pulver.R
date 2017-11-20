#' Computes p-values for interaction terms from linear regressions
#'
#' Given matrices \code{ymat}, \code{xmat}, and \code{zmat},
#' \code{pulverize} evaluates every linear regression
#' \deqn{y = \beta_0 + \beta_1 x + \beta_2 z + \beta_3 xz + \epsilon}{y = b0 + b1 x + b2 z + b3 xz}
#' where \eqn{y}, \eqn{x}, and \eqn{z} are columns of \code{ymat},
#' \code{xmat}, and \code{zmat}, respectively, and returns the p-value
#' for the null hypothesis \eqn{\beta_3=0}{b3=0}.  For example, if
#' \code{ymat}, \code{xmat}, and \code{zmat} have 200, 300, and 400
#' columns, then \code{pulverize} evaluates 24 million models
#' (\eqn{200\times 300\times 400}{200 * 300 * 400}) and returns the
#' p-value for the interaction term for each model.
#'
#' For reasons of computational efficiency, \code{pulverize} returns
#' only p-values and only for the interaction term \eqn{xz}.
#' Fast run time is achieved by using the correlation coefficient between
#' the outcome variable \eqn{y} and the interaction term \eqn{xz} to test
#' the null-hypothesis, which avoids the costly computation of inversions.
#' Additional employed time-saving operations are a rearrangement of the order when iterating through
#' the different matrices, and implementing the core algorithm in the fast
#' programming language C++.
#'
#' Once interesting models are identified based on the resulting p-values, the number of
#' models will be greatly reduced.  At this point additional model
#' characteristics, e.g. effect estimates and standard errors, can be
#' obtained via traditional methods such as R's
#' \code{\link[stats]{lm}} function.
#'
#' @param ymat,xmat,zmat Matrices of type "double", number of rows for all
#'     matrices are the number of observations and must have the same size
#' @param output_file Output file if \code{NULL} no output file will be created
#'     (default: \code{NULL})
#' @param colnames Column names for result table (default: "y", "x", and "z")
#' @param pvalue_threshold Report only p-values below threshold
#'     (default:  \code{NULL}, i.e., report all p-values)
#' @param cores Number of cores to use for parallelization (default:
#'     1)
#' @param overwrite If \code{TRUE} overwrite \code{output_file} (default:
#'     \code{FALSE})
#' @param suppress_return If \code{TRUE} return \code{NULL} instead of a data
#'     frame with p-values for the interaction term (default: \code{FALSE})
#' @param .fp File pointer, only for internal use. Necessary for the function
#'     \code{\link{pulverize_all}} to join results from different input files
#'     to the same output file (dafault: \code{NULL})
#'
#' @details Matrices \code{ymat}, \code{xmat}, and \code{zmat} must
#'     have column names. Missing values are imputed using their
#'     column means. If \code{pvalue_threshold} is supplied, only
#'     p-values (p < \code{pvalue_threshold}) strictly below the
#'     threshold are included in the returned data frame and saved
#'     in \code{output_file}.
#'
#'     In cases where the resulting data frame would be too large to fit
#'     in memory, it is possible to write the results to
#'     \code{output_file} without returning a data frame by setting
#'     \code{suppress_return} to \code{TRUE}.
#'
#'     The column names of the results table are by default set to "y", "x", and
#'     "z", but can be changed using the \code{colnames} argument.
#'
#'     An error will be signaled if \code{output_file} already exists.
#'     Setting \code{overwrite} to \code{TRUE} will silently overwrite the
#'     file.
#'
#'     The computation can be parallelized by specifying a number of
#'     \code{cores} greater than 1.  By default only a single CPU is
#'     used.  Note that parallelization is only supported in
#'     environments with C/C++ compilers that support OpenMP.
#'
#' @return If \code{suppress_return} is \code{FALSE}, a data frame with
#'     columns "y", "x", "z", and "pvalue" containing the p-values for
#'     the interaction term \eqn{xz} of the above linear model is returned.
#'     Otherwise \code{NULL} is returned.
#'
#' @references Shabalin, Andrey A (2012) Bioinformatics: Matrix eQTL:
#'     ultra fast eQTL analysis via large matrix operations, Oxford
#'     Univ Press, *28*, 1353-1358
#'
#' @examples
#' nobs <- 100
#' y <- matrix(rnorm(nobs * 2), ncol = 2, dimnames = list(paste0("row", 1:nobs),
#' paste0("column", 1:2)))
#' x <- matrix(rnorm(nobs * 3), ncol = 3, dimnames = list(paste0("row", 1:nobs),
#' paste0("column", 1:3)))
#' z <- matrix(rnorm(nobs * 4), ncol = 4, dimnames = list(paste0("row", 1:nobs),
#' paste0("column", 1:4)))
#' pulverize(y, x, z)
#'
#' @export
pulverize <- function(ymat, xmat, zmat, output_file = NULL,
    colnames = c("y", "x", "z"), pvalue_threshold = NULL, cores = 1L,
    overwrite = FALSE, suppress_return = FALSE, .fp = NULL)
{
    if (!is.null(output_file) && file.exists(output_file) && !overwrite)
        stop("File exists already: ", output_file)
    if (any(c(is.null(colnames(ymat)), is.null(colnames(xmat)), is.null(colnames(zmat)))))
        stop("Column names for at least one matrix is missing.")
    check_identical_rownames(ymat, xmat, zmat)
    number_of_individuals <- nrow(xmat)
    if (is.null(pvalue_threshold))
        rvalue_threshold <- Inf
    else
        rvalue_threshold <- p2r(pvalue_threshold, df = number_of_individuals - 4)
    swap_x_and_z <- ncol(zmat) > ncol(xmat)
    size <- ncol(ymat) * ncol(if (swap_x_and_z) zmat else xmat)
    storage <- create_storage(size, nrow(ymat), cores)
    on.exit(delete_storage(storage), add = TRUE)
    if (anyNA(ymat))
        ymat <- replace_NA_by_mean(ymat)
    if (anyNA(xmat))
        xmat <- replace_NA_by_mean(xmat)
    if (anyNA(zmat))
        zmat <- replace_NA_by_mean(zmat)
    center <- function(mat) {
        apply(mat, 2, function(x) x - mean(x))
    }
    ymat <- center(ymat)
    xmat <- center(xmat)
    zmat <- center(zmat)
    opened_output_file <- FALSE
    if (is.null(.fp)) {
        .fp <- methods::new("FilePointer", output_file)
        opened_output_file <- TRUE
    }
    result <- tryCatch({
        if (opened_output_file)
            write_header(.fp, c(colnames, "r"))
        .Call('_pulver_compute_and_save_results', PACKAGE = 'pulver',
            ymat, xmat, zmat, rvalue_threshold, colnames,
            min(cores, parallel::detectCores()), suppress_return,
            swap_x_and_z, storage, .fp@pointer)
    }, finally = {
        if (opened_output_file)
            close_output_file(.fp)
    })
    if (!suppress_return)
        result
}

#' pulverize_all
#'
#' Call \code{\link{pulverize}} for all combinations of input files.
#'
#' @param ys,xs,zs lists of input file names (text or \code{databel} object file)
#'     containing tab-separated tables or matrices of type "double"
#' @param output_file Output file name (default: \code{NULL})
#' @param colnames Column names of output file (default: "y", "x", and "z")
#' @param pvalue_threshold p-value threshold (default: NULL, see \code{\link{pulverize}})
#' @param cores Number of cores (default: 1, see \code{\link{pulverize}})
#' @param overwrite If \code{TRUE} overwrite \code{output_file} (default: \code{FALSE})
#' @param suppress_return If TRUE return NULL instead of a data
#'     frame with p-values for the interaction term (default: \code{TRUE})
#'
#' @details \code{pulverize_all} iterates through lists of file names or list
#'     of matrices and calls \code{pulverize} for each combination. Files must contain
#'     matrices with the same number of rows (= number of observations) and must have
#'     column names.
#'
#'     By default all results are written to \code{output_file} without
#'     returning a data frame (\code{suppress_return} = \code{TRUE}). Column names
#'     of the output file are by default set to "y", "x", and "z", but can be changed
#'     using the \code{colnames} argument. Furthermore, by default, \code{output_file} cannot
#'     be overwritten. For enabling overwrite of an existing \code{output_file}, set
#'     \code{overwrite} to \code{TRUE}.
#'
#' @seealso \code{\link{pulverize}}
#'
#' @examples
#' #generate files
#' nobs <- 100
#' y <- matrix(rnorm(nobs * 2), ncol = 2, dimnames = list(NULL,
#' paste0("column", 1:2)))
#' write.table(y, file = "y.txt", sep = "\t")
#' write_databel(y, "y")
#' x <- matrix(rnorm(nobs * 3), ncol = 3, dimnames = list(NULL,
#' paste0("column", 1:3)))
#' write.table(x, file = "x.txt", sep = "\t")
#' z <- matrix(rnorm(nobs * 4), ncol = 4, dimnames = list(NULL,
#' paste0("column", 1:4)))
#' write.table(z, file = "z.txt", sep = "\t")
#'
#' #run pulverize_all with generated files
#' pulverize_all(as.list("y.txt", "y.fvi"), as.list("x.txt"), as.list("z.txt"),
#' suppress_return = FALSE)
#'
#' #remove generated files
#' file.remove("x.txt", "y.txt", "z.txt", "y.fvi", "y.fvd")
#'
#' @export
pulverize_all <- function(ys, xs, zs, output_file = NULL, colnames = c("y", "x", "z"), pvalue_threshold = 1, cores = 1L, overwrite = FALSE, suppress_return = TRUE) {
    if (!is.null(output_file) && file.exists(output_file) && !overwrite)
        stop("File exists already: ", output_file)
    results <- vector("list", length(ys) * length(xs) * length(zs))
    i <- 1L
    fp <- methods::new("FilePointer", output_file)
    tryCatch({
        write_header(fp, colnames)
        for (y in ys) {
            ymat <- if (is.character(y)) read_data(y) else y
            for (x in xs) {
                xmat <- if (is.character(x)) read_data(x) else x
                for (z in zs) {
                    zmat <- if (is.character(z)) read_data(z) else z
                    result <- pulverize(ymat, xmat, zmat, output_file, colnames = colnames,
                        pvalue_threshold = pvalue_threshold, cores = min(cores, parallel::detectCores()),
                        overwrite = TRUE, suppress_return = suppress_return, .fp = fp)
                    results[[i]] <- result
                    i <- i + 1L
                }
            }
        }
    }, finally = close_output_file(fp))
    if (suppress_return)
        NULL
    else
        .Call("_pulver_combine_results", PACKAGE = "pulver", results)
}

setClass("FilePointer", slots = list(pointer = "externalptr"))
setMethod("initialize", "FilePointer", function(.Object, filename) {
    if (is.null(filename))
        .Object@pointer <- .Call("_pulver_create_fake_pointer")
    else
        .Object@pointer <- .Call("_pulver_open_output_file", filename)
    .Object
})

setGeneric("close_output_file", function(fp) stop("not implemented"))
setMethod("close_output_file", "FilePointer", function(fp) {
    .Call("_pulver_close_output_file", fp@pointer)
})

setGeneric("write_header", function(fp, names) stop("not implemented"))
setMethod("write_header", "FilePointer", function(fp, names) {
    .Call("_pulver_write_header", fp@pointer, names)
})

p2t <- function(p, df) {
    stats::qt(p / 2, df, lower.tail = FALSE)
}

t2r <- function(t, df) {
    t * sqrt(1 / (df + t * t))
}

p2r <- function(p, df) {
    t2r(p2t(p, df), df)
}

r2t <- function(r, df) {
    r * sqrt(df / (1 - r * r))
}

read_data <- function(file) {
    if (tools::file_ext(file) == "fvi")
        read_databel(file)
    else
        as.matrix(utils::read.delim(file))
}

#' Reads a \code{databel} object file in table format and returns a matrix of type "double"
#'
#' @param file file name of a \code{databel} object (no extension, i.e. "example" instead of  "example.fvi")
#'
#' @return matrix of type "double"
#'
#' @references \url{www.genabel.org/packages/DatABEL}
#'
#' @export
read_databel <- function(file) {
    x <- DatABEL::databel(tools::file_path_sans_ext(file))
    data <- DatABEL::databel2matrix(x)
    DatABEL::disconnect(x)
    data
}

#' Converts matrix of type "double" to a \code{databel} object file
#'
#' @param data matrix of type "double"
#' @param file file name of a \code{databel} object or text file
#'
#' @references \url{www.genabel.org/packages/DatABEL}
#'
#' @export
write_databel <- function(data, file) {
    x <- DatABEL::matrix2databel(data, file)
    DatABEL::disconnect(x)
    NULL
}

create_storage <- function(size, nrow, cores) {
    .Call("_pulver_create_storage", size, nrow, cores, PACKAGE = "pulver")
}

delete_storage <- function(storage) {
    .Call("_pulver_delete_storage", storage, PACKAGE = "pulver")
}

replace_NA_by_mean <- function(data) {
    datacM <- colMeans(data, na.rm=TRUE)
    indx <- which(is.na(data), arr.ind=TRUE)
    data[indx] <- datacM[indx[,2]]
    data
}

check_identical_rownames <- function(ymat, xmat, zmat) {
    nullmatrices <- c(is.null(rownames(ymat)), is.null(rownames(xmat)), is.null(rownames(zmat)))
    if (any(nullmatrices)) {
        if (length(which(nullmatrices)) == 1) {
            warning(paste0("Matrix ", c("ymat", "xmat", "zmat")[nullmatrices], " does not have row names."))
        } else {
            warning(paste0("Matrices ", paste(c("ymat", "xmat", "zmat")[nullmatrices], collapse = ", ")," do not have row names."))
        }
    } else {
        if (any(rownames(ymat) != rownames(xmat) | rownames(xmat) != rownames(zmat)))
            stop("Matrices do not have identical row names.")
    }
}
