#include "partial_result.h"
#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif
#include <algorithm>
#include <cmath>
#include <numeric>

inline static double r2p(double rvalue, int df);
inline static void center(double *y, double *result, int len);
inline static void ortho(double *y, double *x, double xx, int len);

PartialResult::PartialResult(Storage* storage, int cores, double rvalue_threshold, int swap_x_and_z)
    : storage_(storage), cores_(cores), rvalue_threshold_(rvalue_threshold), swap_x_and_z_(swap_x_and_z) {}

void PartialResult::compute(NumericMatrix ymat, NumericMatrix xmat, NumericMatrix zmat, int zcolumn) {
    run_regressions(ymat, xmat, zmat, zcolumn);
    filter_and_annotate_results(ymat, xmat, zmat, zcolumn);
}

void PartialResult::copy(CharacterVector y, CharacterVector x, CharacterVector z, NumericVector pvalues, int offset) {
    for (int i = 0; i < size(); ++i) {
        y[offset + i] = y_[i];
        x[offset + i] = x_[i];
        z[offset + i] = z_[i];
        pvalues[offset + i] = pvalues_[i];
    }
}

void PartialResult::write(std::FILE *fp) {
    for (int i = 0; i < (int) pvalues_.size(); ++i)
        std::fprintf(fp, "%s\t%s\t%s\t%.17g\n", y_[i], x_[i], z_[i], pvalues_[i]);
}

void PartialResult::run_regressions(NumericMatrix ymat, NumericMatrix xmat, NumericMatrix zmat, int zcolumn) {
    std::vector<double>& rvalues = storage_->rvalues_;
    int n = ymat.nrow();
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores_)
#endif
    for (int j = 0; j < xmat.ncol(); ++j) {
#ifdef _OPENMP
        int thread = omp_get_thread_num();
#else
        int thread = 0;
#endif
        double *y = storage_->y_[thread];
        double *xz = storage_->xz_[thread];
        double *z = storage_->z_[thread];
        double *x = &xmat[j * n];
        double xx = std::inner_product(x, x + n, x, 0.0);
        std::copy(&zmat[zcolumn * n], &zmat[zcolumn * n + n], z);
        for (int m = 0; m < n; ++m)
            xz[m] = x[m] * z[m];
        center(xz, xz, n);
        ortho(z, x, xx, n);
        double zz = std::inner_product(z, z + n, z, 0.0);
        ortho(xz, x, xx, n);
        ortho(xz, z, zz, n);
        double norm_xz = std::sqrt(std::inner_product(xz, xz + n, xz, 0.0));
        for (int k = 0; k < ymat.ncol(); ++k) {
            std::copy(&ymat[k * n], &ymat[k * n] + n, y);
            ortho(y, x, xx, n);
            ortho(y, z, zz, n);
            double norm_y = std::sqrt(std::inner_product(y, y + n, y, 0.0));
            double r = std::inner_product(xz, xz + n, y, 0.0) / (norm_y * norm_xz);
            rvalues[j * ymat.ncol() + k] = r;
        }
    }
}

void PartialResult::filter_and_annotate_results(NumericMatrix ymat, NumericMatrix xmat, NumericMatrix zmat, int zcolumn) {
    std::vector<double>& rvalues = storage_->rvalues_;
    CharacterVector ynames = colnames(ymat);
    CharacterVector xnames = colnames(xmat);
    CharacterVector znames = colnames(zmat);
    int df = ymat.nrow() - 4;
    int no_threshold = !R_FINITE(rvalue_threshold_);
    for (int i = 0; i < (int) rvalues.size(); ++i) {
        if (no_threshold || std::abs(rvalues[i]) > rvalue_threshold_) {
            y_.push_back(CHAR(ynames[i % ymat.ncol()]));
            x_.push_back(CHAR(swap_x_and_z_ ? znames[zcolumn] : xnames[i / ymat.ncol()]));
            z_.push_back(CHAR(swap_x_and_z_ ? xnames[i / ymat.ncol()] : znames[zcolumn]));
            pvalues_.push_back(r2p(rvalues[i], df));
        }
    }
}

inline static double r2p(double r, int df) {
    double t = r * std::sqrt(df / (1 - r*r));
    return 2 * R::pt(std::abs(t), df, 0, 0);
}

inline static void center(double *x, double *result, int len) {
    double mean = std::accumulate(x, x + len, 0.0) / len;
    for (int i = 0; i < len; ++i)
        result[i] = x[i] - mean;
}

inline static void ortho(double *x, double *y, double yy, int len) {
    // Compute x - <x,y> / <y,y> * y.
    double xy = std::inner_product(x, x + len, y, 0.0);
    for (int i = 0; i < len; ++i)
        x[i] -= xy / yy * y[i];
}
