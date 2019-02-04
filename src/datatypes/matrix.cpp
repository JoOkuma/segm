
#include "matrix.h"

using namespace segm;

template<>
void Matrix<float>::mult(const Matrix<float> &B, Matrix<float> &out, bool A_transpose,
                         bool B_transpose, float alpha) const
{
    int A_row = (A_transpose) ? col : row;
    int A_col = (A_transpose) ? row : col;

    int B_row = (B_transpose) ? B.getCol() : B.getRow();
    int B_col = (B_transpose) ? B.getRow() : B.getCol();

    float beta = 0.0;

    if (A_row != out.getRow() || B_col != out.getCol())
        throw std::runtime_error("Dimensions of output matrix (" + std::to_string(out.getRow()) + ", " +
                                 std::to_string(out.getCol()) + ") doesn't match A (" + std::to_string(row) +
                                 ", " + std::to_string(col) + ((A_transpose) ? ")^T" : ")") + "and B (" +
                                 std::to_string(B.getRow()) + ",  " + std::to_string(B.getCol()) +
                                 ((B_transpose) ? ")^T" : ")"));

    if (A_col != B_row)
        throw std::runtime_error("Cannot multiply matrices A is (" + std::to_string(row) + ", " +
                                 std::to_string(col) + ((A_transpose) ? ")^T" : ")") + "and B is (" +
                                 std::to_string(B.getRow()) + ",  " + std::to_string(B.getCol()) +
                                 ((B_transpose) ? ")^T" : ")"));

    cblas_sgemm(CblasRowMajor, (A_transpose) ? CblasTrans : CblasNoTrans,
                (B_transpose)  ? CblasTrans : CblasNoTrans, A_row, B_col,
                A_col, alpha, feat, A_col, B.getFeats(), B_col, beta, out.getFeats(), B_col);
}

template<>
void Matrix<double>::mult(const Matrix<double> &B, Matrix<double> &out, bool A_transpose,
                          bool B_transpose, float alpha) const
{
    int A_row = (A_transpose) ? col : row;
    int A_col = (A_transpose) ? row : col;

    int B_row = (B_transpose) ? B.getCol() : B.getRow();
    int B_col = (B_transpose) ? B.getRow() : B.getCol();

    float beta = 0.0;

    if (A_row != out.getRow() || B_col != out.getCol())
        throw std::runtime_error("Dimensions of output matrix (" + std::to_string(out.getRow()) + ", " +
                                 std::to_string(out.getCol()) + ") doesn't match A (" + std::to_string(row) +
                                 ", " + std::to_string(col) + ((A_transpose) ? ")^T" : ")") + "and B (" +
                                 std::to_string(B.getRow()) + ",  " + std::to_string(B.getCol()) +
                                 ((B_transpose) ? ")^T" : ")"));

    if (A_col != B_row)
        throw std::runtime_error("Cannot multiply matrices A is (" + std::to_string(row) + ", " +
                                 std::to_string(col) + ((A_transpose) ? ")^T" : ")") + "and B is (" +
                                 std::to_string(B.getRow()) + ",  " + std::to_string(B.getCol()) +
                                 ((B_transpose) ? ")^T" : ")"));

    cblas_dgemm(CblasRowMajor, (A_transpose) ? CblasTrans : CblasNoTrans,
                (B_transpose)  ? CblasTrans : CblasNoTrans, A_row, B_col,
                A_col, alpha, feat, A_col, B.getFeats(), B_col, beta, out.getFeats(), B_col);
}
