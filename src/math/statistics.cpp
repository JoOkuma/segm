
#include "statistics.h"

using namespace Eigen;

namespace segm {

    template <typename T>
    Matrix<T, Dynamic, Dynamic> PCA(Matrix<T, Dynamic, Dynamic> &data, int dimension_out)
    {
        long dim_in = data.cols();
        long n = data.rows();

        if (dimension_out > dim_in)
            throw std::invalid_argument("PCA output dimension cannot be greater than input dimension");

        Matrix<T, Dynamic, 1> mean(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < dim_in; j++) {
                mean[j] += data(i, j);
            }
        }

        mean /= n;

        Matrix<T, Dynamic, Dynamic> centralized(n, dim_in);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < dim_in; j++) {
                centralized(i, j) = data(i, j) - mean[j];
            }
        }

//        Matrix<T, Dynamic, Dynamic> cov = centralized.transpose() * centralized();
        Matrix<T, Dynamic, Dynamic> cov = centralized.transpose().lazyProduct(centralized);
        cov *= (1.0 / (n - 1));

        BDCSVD<Matrix<T, Dynamic, Dynamic>> svd(cov, ComputeFullV);

        Matrix<T, Dynamic, Dynamic> V = svd.matrixV();
        V.conservativeResize(dimension_out, dim_in);

        return V;
    }

    VectorXi randomIndexes(int size) {
        VectorXi index(size);
        for (int i = 0; i < size; i++)
            index[i] = i;

        for (int i = 0; i < size; i++) {
            int j = (int) (i + random() / (RAND_MAX / (size - 1) + 1));
            int tmp = index[j];
            index[j] = index[i];
            index[i] = tmp;
        }
        return index;
    }


    template Matrix<float, Dynamic, Dynamic> PCA(Matrix<float, Dynamic, Dynamic> &data, int dimension_out);
    template Matrix<double, Dynamic, Dynamic> PCA(Matrix<double, Dynamic, Dynamic> &data, int dimension_out);
}