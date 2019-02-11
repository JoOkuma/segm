
#include "math/statistics.h"
#include "pcatest.h"
#include "test.h"

#include <cmath>

void pca_test()
{
    Eigen::MatrixXd input(4, 2);
    input << 5, 2,
             3, 4,
             5, 6,
             7, 8;

    Eigen::MatrixXd pca = segm::PCA(input, 1);

    double res_pca[2] = {0.4472136, 0.89442719};

    for (int i = 0; i < pca.cols(); i++) {
        ASSERT_THROW((fabs(pca(0, i) - res_pca[i]) < 1e-5))
    }
}