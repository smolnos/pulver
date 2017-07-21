library(methods)
library(pulver)
source("pulverize_mixabel.R")
source("pulverize_matrixeqtl.R")

set_parameters <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    params <- setNames(as.list(argv[c(FALSE, TRUE)]), argv[c(TRUE, FALSE)])
    argnames <- c("nobs", "xcols", "ycols", "zcols", "cores", "code", "tag")
    if(any(!argnames %in% names(params)))
        stop("Usage: Rscript run_benchmark.R nobs N xcols N ycols N zcols N cores N code \"function\" tag \"description\"")
    numbers <- setdiff(argnames, c("code", "tag"))
    params[numbers] <- lapply(params[numbers], as.integer)
    params["commit_id"] <- substr(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), 1, 6)
    params["seed"] <- 1234
    params
}

summarize_benchmark <- function(timing, params) {
    data.frame(
        tag = params$tag,
        commit_id = params$commit_id,
        ddmmyy = format(Sys.time(), "%d/%m/%y"),
        time =  format(Sys.time(), "%H:%M"),
        code = params$code,
        nobs = params$nobs,
        xcols = params$xcols,
        ycols = params$ycols,
        zcols = params$zcols,
        cores = params$cores,
        mean_runtime = mean(timing$time) / 1e6,
        sd_runtime = sd(timing$time) / 1e6,
        runs = nrow(timing),
        stringsAsFactors = FALSE)
}

params <- set_parameters()

set.seed(params$seed)

mat <- lapply(pulver:::create_input_files(
    dir = "../data",
    type = "fvi",
    nobs = params$nobs,
    xcols = params$xcols,
    ycols = params$ycols,
    zcols = params$zcols), pulver:::read_data)

output_file <- sprintf("benchmark_%s_seed%d_function%s_nobs%d_xcols%d_ycols%d_zcols%d_cores%d.txt",
    params$commit_id, params$seed, params$code, params$nobs, params$xcols, params$ycols, params$zcols, params$cores)
log_file <- paste0(tools::file_path_sans_ext(output_file), ".log")

code <- switch(params$code,
    pulverize = expression(pulverize(mat$x, mat$y, mat$z, output_file, cores = params$cores, overwrite = TRUE)),
    mixabel = expression(pulverize_mixabel(mat$x, mat$y, mat$z, pvalue_threshold = 1, colnames = c("x", "y", "z"))),
    matrixeqtl = expression(pulverize_matrixeqtl(mat$x, mat$y, mat$z, output_file, pvalue_threshold = 1)),
    lm = expression(pulver:::pulverize_lm(mat$x, mat$y, mat$z, pvalue_threshold = 1, colnames = c("x", "y", "z"))),
    foydle = expression(foydle::foydle(mat$x, mat$y, mat$z, output_file, cores = params$cores))
)

print(params$code)
timing <- microbenchmark::microbenchmark(eval.parent(code), times = 200)

summary <- summarize_benchmark(timing, params)
write.table(summary, log_file, quote = FALSE, sep = "\t", row.names = FALSE)
print(summary[c("mean_runtime", "sd_runtime")])
