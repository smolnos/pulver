#include "../inst/include/pulver.h"

Storage::Storage(int size, int nrow, int cores) : cores_(cores) {
    rvalues_.resize(size);
    for (int i = 0; i < cores_; ++i) {
        y_.push_back(new double[nrow]);
        z_.push_back(new double[nrow]);
        xz_.push_back(new double[nrow]);
    }
}

Storage::~Storage() {
    for (int i = 0; i < cores_; ++i) {
        delete [] y_[i];
        delete [] z_[i];
        delete [] xz_[i];
    }
}
