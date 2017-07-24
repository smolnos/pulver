#include <RcppCommon.h>
#include "../inst/include/pulver.h"
#include <Rcpp.h>
#include "partial_result.h"
#include "result.h"
#include <cstdio>

#ifdef SUPPORT_OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;

static int count_rows(List results);

// [[Rcpp::export]]
DataFrame compute_and_save_results(NumericMatrix ymat, SEXP xmat_,
    SEXP zmat_, double rvalue_threshold, CharacterVector names,
    int cores, int suppress_return, int swap_x_and_z, SEXP xpstorage, SEXP xp)
{
    NumericMatrix xmat(swap_x_and_z ? zmat_ : xmat_);
    NumericMatrix zmat(swap_x_and_z ? xmat_ : zmat_);
    FILE *fp = (FILE *) R_ExternalPtrAddr(xp);
    XPtr<Storage> ptr(xpstorage);
    Storage *storage = ptr;
    Result result;
    for (int zcolumn = 0; zcolumn < zmat.ncol(); ++zcolumn) {
        PartialResult *partial_result = new PartialResult(storage, cores, rvalue_threshold, swap_x_and_z);
        partial_result->compute(ymat, xmat, zmat, zcolumn);
        checkUserInterrupt();
        if (fp)
            partial_result->write(fp);
        if (suppress_return)
            delete partial_result;
        else
            result.add(partial_result);
    }
    if (suppress_return)
        return R_NilValue;
    return result.as_data_frame(names);
}

// [[Rcpp::export]]
void write_header(SEXP xp, CharacterVector names) {
    FILE *fp = (FILE *) R_ExternalPtrAddr(xp);
    if (fp)
        std::fprintf(fp, "%s\t%s\t%s\tpvalue\n", CHAR(names[0]), CHAR(names[1]), CHAR(names[2]));
}

// [[Rcpp::export]]
SEXP open_output_file(CharacterVector filename)
{
    const char *file = CHAR(filename[0]);
    FILE *fp = std::fopen(file, "wb");
    return R_MakeExternalPtr(fp, R_NilValue, R_NilValue);
}

// [[Rcpp::export]]
void close_output_file(SEXP xp)
{
    FILE *fp = (FILE *) R_ExternalPtrAddr(xp);
    if (fp)
        std::fclose(fp);
}

// [[Rcpp::export]]
SEXP create_fake_pointer() {
    return R_MakeExternalPtr(0, R_NilValue, R_NilValue);
}

// [[Rcpp::export]]
DataFrame combine_results(List results) {
    int nrow = count_rows(results);
    CharacterVector y(nrow);
    CharacterVector x(nrow);
    CharacterVector z(nrow);
    NumericVector pvalues(nrow);
    int row = 0;//sdf
    for (int i = 0; i < results.size(); ++i) {
        DataFrame d = static_cast<DataFrame> (results[i]);
        CharacterVector y_ = d[0];
        CharacterVector x_ = d[1];
        CharacterVector z_ = d[2];
        NumericVector pvalues_ = d[3];
        for (int j = 0; j < y_.size(); ++j) {
            y[row] = y_[j];
            x[row] = x_[j];
            z[row] = z_[j];
            pvalues[row] = pvalues_[j];
            ++row;
        }
    }
    DataFrame d = static_cast<DataFrame> (results[0]);
    CharacterVector names = d.attr("names");
    return DataFrame::create(
        Named(CHAR(names[0])) = y,
        Named(CHAR(names[1])) = x,
        Named(CHAR(names[2])) = z,
        Named(CHAR(names[3])) = pvalues);
}

static int count_rows(List results) {
    int n = 0;
    for (int i = 0; i < results.size(); ++i) {
        DataFrame d = static_cast<DataFrame> (results[i]);
        CharacterVector y = d[0];
        n += y.size();
    }
    return n;
}

// [[Rcpp::export]]
XPtr<Storage> create_storage(int size, int nrow, int cores){
    Storage *storage = new Storage(size, nrow, cores);
    XPtr<Storage> ptr(storage, false);
    return ptr;
}

// [[Rcpp::export]]
void delete_storage(SEXP xp) {
    XPtr<Storage> ptr(xp);
    Storage *storage = ptr;
    delete storage;
}
