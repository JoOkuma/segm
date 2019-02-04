
#ifndef SEGM_STATISTICS_H
#define SEGM_STATISTICS_H

#include "math/linearalgebra.h"

namespace segm {

    /**
     * @param data              Input matrix, (m, n)
     * @param dimension_out     Number of output dimensions
     * @return                  PCA rotation matrix, (dimension_out, n) with principal components vectors on rows
     */

    Matrix<double> PCA(Matrix<double> &data, int dimension_out);

    Matrix<float> PCA(Matrix<float> &data, int dimension_out);

    Vector<int> randomIndexes(int size);

}

#endif //SEGM_STATISTICS_H
