
#include "linearalgebra.h"

extern "C"
{
extern void dgesdd_(char *jobz, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt,
                    int *ldvt, double *work, int *lwork, int *iwork, int *info);

extern void sgesdd_(char *jobz, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt,
                    int *ldvt, float *work, int *lwork, int *iwork, int *info);
}

namespace segm {


    void SVD(const Matrix<double> &matrix, Matrix<double> &U, Vector<double> &S, Matrix<double> &Vt)
    {
        // A = U*S*Vt
        int A_row = matrix.getRow();
        int A_col = matrix.getCol();
        int A_min_dim = ((A_row < A_col) ? A_row : A_col);

        if (A_row != U.getRow() || A_row != U.getCol())
            throw std::invalid_argument("Matrix U must be a squared matrix with same dimensions"
                                        "as the number of rows of the input");


        if (A_col != Vt.getRow() || A_col != Vt.getCol())
            throw std::invalid_argument("Matrix Vt must be a squared matrix with same dimensions"
                                        "as the number of columns of the input");

        if (S.getSize() != A_min_dim)
            throw std::invalid_argument("Vector S must have the same length as the smaller dimension of the input");

        char jobz = 'A';

        int lwork = -1;
        double work_param;
        int *iwork = new int[A_min_dim * 8];
        int info;

        // lapack is column major
        Matrix<double> A_t = matrix.t();

        dgesdd_(&jobz, &A_row, &A_col, A_t.getFeats(), &A_row, S.getFeats(), U.getFeats(), &A_row,
                Vt.getFeats(), &A_col, &work_param, &lwork, iwork, &info);

        lwork = (int) work_param;
        double *work = new double[lwork];

        dgesdd_(&jobz, &A_row, &A_col, A_t.getFeats(), &A_row, S.getFeats(), U.getFeats(), &A_row,
                Vt.getFeats(), &A_col, work, &lwork, iwork, &info);

        if (info != 0)
            throw std::runtime_error("Eigen matrix computation was not successful");

        delete[] iwork;
        delete[] work;

        Vt = Vt.t();
        U = U.t();
    }


    void SVD(const Matrix<float> &matrix, Matrix<float> &U, Vector<float> &S, Matrix<float> &Vt)
    {
        // A = U*S*Vt
        int A_row = matrix.getRow();
        int A_col = matrix.getCol();
        int A_min_dim = ((A_row < A_col) ? A_row : A_col);

        if (A_row != U.getRow() || A_row != U.getCol())
            throw std::invalid_argument("Matrix U must be a squared matrix with same dimensions"
                                        "as the number of rows of the input");


        if (A_col != Vt.getRow() || A_col != Vt.getCol())
            throw std::invalid_argument("Matrix Vt must be a squared matrix with same dimensions"
                                        "as the number of columns of the input");

        if (S.getSize() != A_min_dim)
            throw std::invalid_argument("Vector S must have the same length as the smaller dimension of the input");

        char jobz = 'A';
        int lwork = -1;
        float work_param;
        int *iwork = new int[A_min_dim * 8];
        int info;

        // lapack is column major
        Matrix<float> A_t = matrix.t();

        sgesdd_(&jobz, &A_row, &A_col, A_t.getFeats(), &A_row, S.getFeats(), U.getFeats(), &A_row,
                Vt.getFeats(), &A_col, &work_param, &lwork, iwork, &info);

        lwork = (int) work_param;
        float *work = new float[lwork];

        sgesdd_(&jobz, &A_row, &A_col, A_t.getFeats(), &A_row, S.getFeats(), U.getFeats(), &A_row,
                Vt.getFeats(), &A_col, work, &lwork, iwork, &info);

        if (info != 0)
            throw std::runtime_error("Eigen matrix computation was not successful");

        delete[] iwork;
        delete[] work;

        Vt = Vt.t();
        U = U.t();
    }

}
