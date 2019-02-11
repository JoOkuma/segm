
#ifndef SEGM_STATISTICS_H
#define SEGM_STATISTICS_H

#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace segm {

    /**
     * @param data              Input matrix, (m, n)
     * @param dimension_out     Number of output dimensions
     * @return                  PCA rotation matrix, (dimension_out, n) with principal components vectors on rows
     */

    template <typename T>
    Matrix<T, Dynamic, Dynamic> PCA(Matrix<T, Dynamic, Dynamic> &data, int dimension_out);

    VectorXi randomIndexes(int size);

}

#endif //SEGM_STATISTICS_H
