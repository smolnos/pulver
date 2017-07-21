#ifndef PULVER_H
#define PULVER_H

#include <vector>

class Storage {
  public:
    Storage(int size, int nrow, int cores);
    ~Storage();
    std::vector<double> rvalues_;
    std::vector<double*> y_;
    std::vector<double*> z_;
    std::vector<double*> xz_;
    int cores_;
};

#endif
