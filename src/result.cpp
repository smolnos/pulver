#include "result.h"

Result::~Result() {
    for (int i = 0; i < (int) results_.size(); ++i)
        delete results_[i];
}

void Result::add(PartialResult* partial_result) {
    results_.push_back(partial_result);
}

DataFrame Result::as_data_frame(CharacterVector names) {
    int nrow = size();
    CharacterVector y(nrow);
    CharacterVector x(nrow);
    CharacterVector z(nrow);
    NumericVector pvalues(nrow);
    int copied = 0;
    for (int i = 0; i < (int) results_.size(); ++i) {
        results_[i]->copy(y, x, z, pvalues, copied);
        copied += results_[i]->size();
    }
    return DataFrame::create(
        Named(CHAR(names[0])) = y,
        Named(CHAR(names[1])) = x,
        Named(CHAR(names[2])) = z,
        Named("pvalue") = pvalues);
}

int Result::size() {
    int n = 0;
    for (int i = 0; i < (int) results_.size(); ++i)
        n += results_[i]->size();
    return n;
}
