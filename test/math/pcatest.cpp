
#include "math/statistics.h"
#include "pcatest.h"
#include "test.h"

#include <cmath>

void pca_test()
{
    double array[8] = {5, 2,
                       3, 4,
                       5, 6,
                       7, 8};

    segm::Matrix<double> input(4, 2, array);

    segm::Matrix<double> pca = segm::PCA(input, 1);

    double array_pca[2] = {-0.4472136, -0.89442719};

    segm::Matrix<double> res_pca(1, 2, array_pca);

    for (int i = 0; i < pca.getRow(); i++) {
        for (int j = 0; j < pca.getRow(); j++) {
            ASSERT_THROW((fabs(pca(i, j) - res_pca(i, j)) < 1e-5))
        }
    }
}