create_matrix <- function(nrow, ncol, column_prefix, column_offset = 0L) {
    matrix(stats::rnorm(nrow * ncol), nrow, ncol,
        dimnames = list(1:nrow, paste0(column_prefix, column_offset + 1:ncol)))
}

create_multiple_xyz_matrices <- function(y = 1, x = 1, z = 1, ycols = 2, xcols = 3, zcols = 4, nobs = 10) {
    ys <- lapply(1:y, function(k) create_matrix(nobs, ycols, "y", (k-1) * ycols))
    xs <- lapply(1:x, function(k) create_matrix(nobs, xcols, "x", (k-1) * xcols))
    zs <- lapply(1:z, function(k) create_matrix(nobs, zcols, "z", (k-1) * zcols))
    list(y = ys, x = xs, z = zs)
}

create_xyz_matrices <- function(ycols = 2, xcols = 3, zcols = 4, nobs = 10) {
    lapply(create_multiple_xyz_matrices(ycols = ycols, xcols = xcols, zcols = zcols, nobs = nobs), `[[`, 1)
}

create_xyz_files <- function(y = 1, x = 1, z = 1, ycols = 2, xcols = 3, zcols = 4,
    type = "txt", dir = "data", nobs = 10) {
    type_index <- 1L
    create_filenames <- function(n) {
        root <- deparse(substitute(n))
        sapply(1:n, function(k) {
            result <- file.path(dir, paste0(root, k, ".", type[type_index]))
            type_index <<- (type_index %% length(type)) + 1
            result
        })
    }
    files <- list(
        y = create_filenames(y),
        x = create_filenames(x),
        z = create_filenames(z))
    Map(write_datasets, files, list(ycols, xcols, zcols), list("y", "x", "z"), nobs)
    files
}

write_datasets <- function(files, ncol = 2, column_prefix = "col", nobs = 10) {
    column <- 1L
    lapply(files, function(file) {
        column_suffix <- seq(column, column + ncol - 1L)
        column <<- column + ncol
        data <- create_matrix(nobs, ncol, column_prefix, column)
        unlink(file, force = TRUE)
        if (tools::file_ext(file) == "fvi")
            write_databel(data, tools::file_path_sans_ext(file))
        else
            utils::write.table(data, file, sep = "\t")
    })
    NULL
}

pulverize_lm_all <- function(ys, xs, zs, pvalue_threshold = 1, colnames = c("y", "x", "z")) {
    result <- vector("list", length(ys) * length(xs) * length(zs))
    i <- 1L
    for (y in ys) {
        ymat <- if (is.character(y)) read_data(y) else y
        for (x in xs) {
            xmat <- if (is.character(x)) read_data(x) else x
            for (z in zs) {
                zmat <- if (is.character(z)) read_data(z) else z
                result[[i]] <- pulverize_lm(ymat, xmat, zmat, pvalue_threshold, colnames)
                i <- i + 1L
            }
        }
    }
    do.call(rbind, result)
}

pulverize_lm <- function(ymat, xmat, zmat, pvalue_threshold = 1, colnames = c("y", "x", "z")) {
    result <- cbind(create_pulverize_annotation(ymat, xmat, zmat, colnames), pvalue = NA_real_)
    for (i in colnames(ymat))
        for (j in colnames(xmat))
            for (k in colnames(zmat)) {
                y <- ymat[, i]
                x <- xmat[, j]
                z <- zmat[, k]
                pvalue <- stats::coef(summary(stats::lm(y ~ x * z)))["x:z", "Pr(>|t|)"]
                result$pvalue[with(result, y == i & x == j & z == k)] <- pvalue
            }
    result[which(result$pvalue < pvalue_threshold), ]
}

create_pulverize_annotation <- function(y, x, z, colnames = c("y", "x", "z")) {
   annot <- data.frame(
        y = rep(colnames(y), each = ncol(x) * ncol(z)),
        x = rep(rep(colnames(x), each = ncol(z)), times = ncol(y)),
        z = rep(colnames(z), times = ncol(y) * ncol(x)),
        stringsAsFactors = FALSE)
   colnames(annot) <- colnames
   annot
}

r2t <- function(r, df) {
    r * sqrt(df / (1 - r * r))
}

r2p <- function(r, df) {
    t2p(r2t(r, df), df)
}

t2p <- function(t, df) {
    2 * stats::pt(abs(t), df, lower.tail = FALSE)
}

p2r <- function(p, df) {
    t2r(p2t(p, df), df)
}

p2t <- function(p, df) {
    stats::qt(p / 2, df, lower.tail = FALSE)
}

t2r <- function(t, df) {
    t * sqrt(1 / (df + t * t))
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

replace_NA_by_mean <- function(data) {
    datacM <- colMeans(data, na.rm=TRUE)
    indx <- which(is.na(data), arr.ind=TRUE)
    data[indx] <- datacM[indx[,2]]
    data
}
