
#ifndef SEGM_ALGORITHM_H
#define SEGM_ALGORITHM_H

#include <vector>
#include <numeric>
#include <algorithm>

template<typename T>
std::vector<int> indexesDecreasing(const std::vector<T> &v) {

    std::vector<int> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

    return idx;
}


template<typename T>
std::vector<int> indexesIncreasing(const std::vector<T> &v) {

    std::vector<int> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}



#endif //SEGM_ALGORITHM_H
