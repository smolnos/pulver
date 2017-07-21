#ifndef RESULT_H
#define RESULT_H

#include "Rcpp.h"
#include "partial_result.h"
#include <vector>

using namespace Rcpp;

class Result {
  public:
    ~Result();
    void add(PartialResult*);
    DataFrame as_data_frame(CharacterVector names);
  private:
    std::vector<PartialResult*> results_;
    int size();
};

#endif
