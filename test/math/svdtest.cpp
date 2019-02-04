
#include "math/linearalgebra.h"
#include "svdtest.h"
#include "test.h"

#include <cmath>

void svd_test()
{
    double array[8] = {1, 2, 3, 4,
                       5, 6, 7, 8};

    segm::Matrix<double> input(2, 4, array);

    segm::Matrix<double> U(input.getRow(), input.getRow());
    segm::Vector<double> S(input.getRow());
    segm::Matrix<double> Vt(input.getCol(), input.getCol());

    segm::SVD(input, U, S, Vt);

    // from numpy
    double array_U[4] = {-0.376168, -0.926551,
                         -0.926551, 0.376168};

    double array_S[2] = {14.2274, 1.25733};

    double array_Vt[16] = {-0.35206169, -0.44362578, -0.53518987, -0.62675396,
                            0.75898127,  0.3212416 , -0.11649807, -0.55423774,
                           -0.40008743,  0.25463292,  0.69099646, -0.54554195,
                           -0.37407225,  0.79697056, -0.47172438,  0.04882607};

    segm::Matrix<double> res_U(input.getRow(), input.getRow(), array_U);
    segm::Vector<double> res_S(input.getRow(), array_S);
    segm::Matrix<double> res_Vt(input.getCol(), input.getCol(), array_Vt);

    for (int i = 0; i < input.getRow(); i++) {
        for (int j = 0; j < input.getRow(); j++) {
            ASSERT_THROW((fabs(U(i, j) - res_U(i, j)) < 1e-5))
        }
    }

    for (int i = 0; i < input.getRow(); i++) {
        ASSERT_THROW((fabs(S[i] - res_S[i]) < 1e-5))
    }

    for (int i = 0; i < input.getCol(); i++) {
        for (int j = 0; j < input.getCol(); j++) {
            ASSERT_THROW((fabs(Vt(i, j) - res_Vt(i, j)) < 1e-5))
        }
    }
}