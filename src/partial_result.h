#ifndef PARTIAL_RESULT_H
#define PARTIAL_RESULT_H

#include <Rcpp.h>
#include <cstdio>
#include <vector>
#include "../inst/include/pulver.h"

using namespace Rcpp;

class PartialResult {
  public:
    PartialResult(Storage*, int cores, double rvalue_threshold, int swap_x_and_z);
    void compute(NumericMatrix ymat, NumericMatrix xmat, NumericMatrix zmat, int zcolumn);
    void copy(CharacterVector y, CharacterVector x, CharacterVector z, NumericVector pvalues, int offset);
    int size() { return pvalues_.size(); }
    int empty() { return size() == 0; }
    void write(std::FILE *fp);
  private:
    Storage* storage_;
    int cores_;
    double rvalue_threshold_;
    int swap_x_and_z_;
    std::vector<double> pvalues_;
    std::vector<const char *> y_, x_, z_;
    void run_regressions(NumericMatrix ymat, NumericMatrix xmat, NumericMatrix zmat, int zcolumn);
    void filter_and_annotate_results(NumericMatrix ymat, NumericMatrix xmat, NumericMatrix zmat, int zcolumn);
};

#endif
