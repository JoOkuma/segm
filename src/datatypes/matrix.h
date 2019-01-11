//
// Created by jookuma on 11/01/19.
//

#ifndef SEGM_MATRIX_H
#define SEGM_MATRIX_H

#include <stdexcept>
#include <string>

#include "cblas.h"

namespace segm
{
    template <typename T>
    class Matrix
    {
    public:
        Matrix(int _row, int _col);
        ~Matrix();

        Matrix<T> mult(Matrix<T> &B, bool A_transpose = false, bool B_transpose = false, float alpha = 1.0);
        Matrix<T> operator*(Matrix<T> &B);


        int getRow() const { return row; }
        int getCol() const { return col; }
        T* getFeat() const { return feat; }

        T &operator()(int _row, int _col) { return feat[row_index[_row] + _col]; }

    protected:
        int row;
        int col;

        T *feat = nullptr;
        int *row_index = nullptr;
    };

    template <typename T>
    Matrix<T>::Matrix(int _row, int _col) : row(_row), col(_col)
    {
        feat = new T[row * col];
        row_index = new int[row];
        for (int i = 0, r = 0; i < row; i++, r += col)
            row_index[i] = r;
    }

    template <typename T>
    Matrix<T>::~Matrix()
    {
        delete feat;
        delete row_index;
    }


    template <>
    Matrix<float> Matrix<float>::mult(Matrix<float> &B, bool A_transpose, bool B_transpose, float alpha)
    {
        int A_row = (A_transpose) ? col : row;
        int A_col = (A_transpose) ? row : col;

        int B_row = (B_transpose) ? col : row;
        int B_col = (B_transpose) ? row : col;

        float beta = 0.0;

        if (A_col != B_row)
            std::runtime_error("Cannot multiply matrices A is (" + std::to_string(row) + ", " +
                               std::to_string(col) + ((A_transpose) ? ")^T" : ")") + "and B is (" +
                               std::to_string(B.getRow()) + ",  " + std::to_string(B.getCol()) +
                               ((B_transpose) ? ")^T" : ")"));

        Matrix<float> C(A_row, B_col);

        cblas_sgemm(CblasRowMajor, (A_transpose) ? CblasTrans : CblasNoTrans,
                    (B_transpose)  ? CblasTrans : CblasNoTrans, A_row, B_col,
                    A_col, alpha, feat, A_col, B.getFeat(), B_col, beta, C.getFeat(), B_col);

        return C;
    }

    template <>
    Matrix<float> Matrix<float>::operator*(Matrix<float> &B) {
        return mult(B);
    }

}


#endif //SEGM_MATRIX_H
