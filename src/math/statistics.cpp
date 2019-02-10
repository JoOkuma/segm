
#include "statistics.h"

namespace segm {

    Matrix<double> PCA(Matrix<double> &data, int dimension_out) {
        int dim_in = data.getCol();
        int n = data.getRow();

        if (dimension_out > dim_in)
            throw std::invalid_argument("PCA output dimension cannot be greater than input dimension");

        Vector<double> mean(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < dim_in; j++) {
                mean[j] += data(i, j);
            }
        }

        mean /= n;

        Matrix<double> centralized(n, dim_in);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < dim_in; j++) {
                centralized(i, j) = data(i, j) - mean[j];
            }
        }

        Matrix<double> cov = centralized.mult(centralized, true, false, 1.0f / (n - 1));

        Vector<double> S(cov.getRow());
        Matrix<double> U(cov.getRow(), cov.getRow());
        Matrix<double> Vt(cov.getCol(), cov.getCol());

        SVD(cov, U, S, Vt);

        return Vt.trimRows(0, dimension_out);
    }


    Matrix<float> PCA(Matrix<float> &data, int dimension_out) {
        int dim_in = data.getCol();
        int n = data.getRow();

        if (dimension_out > dim_in)
            throw std::invalid_argument("PCA output dimension cannot be greater than input dimension");

        Vector<float> mean(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < dim_in; j++) {
                mean[j] += data(i, j);
            }
        }

        mean /= n;

        Matrix<float> centralized(n, dim_in);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < dim_in; j++) {
                centralized(i, j) = data(i, j) - mean[j];
            }
        }

        Matrix<float> cov = centralized.mult(centralized, true, false, 1.0f / (n - 1));

        Vector<float> S(cov.getRow());
        Matrix<float> U(cov.getRow(), cov.getRow());
        Matrix<float> Vt(cov.getCol(), cov.getCol());

        SVD(cov, U, S, Vt);

        return Vt.trimRows(0, dimension_out);
    }


    Vector<int> randomIndexes(int size) {
        Vector<int> index(size);
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

}