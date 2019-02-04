
#ifndef SEGM_LINEARALGEBRA_H
#define SEGM_LINEARALGEBRA_H

#include "datatypes/matrix.h"

namespace segm {

    /**
     * @param matrix        Input matrix, (m, n)
     * @param U             Output matrix U, (m, m)
     * @param S             Output vector S, min(m, n)
     * @param Vt            Output transpose matrix V, (n, n)
     */

    void SVD(const Matrix<double> &matrix, Matrix<double> &U, Vector<double> &S, Matrix<double> &Vt);

    void SVD(const Matrix<float> &matrix, Matrix<float> &U, Vector<float> &S, Matrix<float> &Vt);

}

#endif //SEGM_LINEARALGEBRA_H
