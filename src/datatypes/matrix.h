
#ifndef SEGM_MATRIX_H
#define SEGM_MATRIX_H

#include <iostream>
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
        Matrix(int _row, int _col, const T *array);
        ~Matrix();

        void mult(const Matrix<T> &B, Matrix<T> &out, bool A_transpose = false,
                  bool B_transpose = false, float alpha = 1.0) const
        {
            int A_row = (A_transpose) ? col : row;
            int A_col = (A_transpose) ? row : col;

            int B_row = (B_transpose) ? B.getCol() : B.getRow();
            int B_col = (B_transpose) ? B.getRow() : B.getCol();

            float beta = 0.0;

            if (A_row != out.getRow() || B_col != out.getCol())
                std::runtime_error("Dimensions of output matrix (" + std::to_string(out.getRow()) + ", " +
                                   std::to_string(out.getCol()) + ") doesn't match A (" + std::to_string(row) +
                                   ", " + std::to_string(col) + ((A_transpose) ? ")^T" : ")") + "and B (" +
                                   std::to_string(B.getRow()) + ",  " + std::to_string(B.getCol()) +
                                   ((B_transpose) ? ")^T" : ")"));

            if (A_col != B_row)
                std::runtime_error("Cannot multiply matrices A is (" + std::to_string(row) + ", " +
                                   std::to_string(col) + ((A_transpose) ? ")^T" : ")") + "and B is (" +
                                   std::to_string(B.getRow()) + ",  " + std::to_string(B.getCol()) +
                                   ((B_transpose) ? ")^T" : ")"));

            cblas_sgemm(CblasRowMajor, (A_transpose) ? CblasTrans : CblasNoTrans,
                        (B_transpose)  ? CblasTrans : CblasNoTrans, A_row, B_col,
                        A_col, alpha, feat, A_col, B.getFeat(), B_col, beta, out.getFeat(), B_col);


        }

        Matrix<T> mult(const Matrix<T> &B, bool A_transpose = false,
                       bool B_transpose = false, float alpha = 1.0) const;

        Matrix<T> operator*(Matrix<T> &B);

        int getRow() const { return row; }

        int getCol() const { return col; }

        T* getFeat() const { return feat; }

        T* getFeat(int i) const { return &feat[row_index[i]]; }

        void setFeat(const T *array);

        T &operator()(int _row, int _col) const { return feat[row_index[_row] + _col]; }

        T &operator()(int idx) const { return feat[idx]; }

        void print();

        Matrix<T> copy();

    protected:
        int row;
        int col;

        T *feat = nullptr;
        int *row_index = nullptr;
    };

    template <typename T>
    Matrix<T>::Matrix(int _row, int _col) : row(_row), col(_col)
    {
        feat = new T[row * col]();
        row_index = new int[row];
        for (int i = 0, r = 0; i < row; i++, r += col)
            row_index[i] = r;
    }

    template <typename T>
    Matrix<T>::Matrix(int _row, int _col, const T *array) : row(_row), col(_col)
    {
        feat = new T[row * col]();
        row_index = new int[row];
        for (int i = 0, r = 0; i < row; i++, r += col)
            row_index[i] = r;

        setFeat(array);
    }

    template <typename T>
    Matrix<T>::~Matrix()
    {
        delete[] feat;
        delete[] row_index;
    }

    template <typename T>
    void Matrix<T>::setFeat(const T *array) {
        for (int i = 0; i < row * col; i++)
            feat[i] = array[i];
    }

    template <typename T>
    Matrix<T> Matrix<T>::mult(const Matrix<T> &B, bool A_transpose,
                                      bool B_transpose, float alpha) const
    {
        int A_row = (A_transpose) ? col : row;
        int B_col = (B_transpose) ? B.getRow() : B.getCol();

        Matrix<float> C(A_row, B_col);

        mult(B, C, A_transpose, B_transpose, alpha);

        return C;
    }

    template <typename T>
    Matrix<T> Matrix<T>::operator*(Matrix<T> &B) {
        return mult(B);
    }

    template <typename T>
    void Matrix<T>::print() {
        for (int i = 0; i < row; i++) {
            std::cout << (*this)(i, 0);
            for (int j = 1; j < col; j++) {
                std::cout << " " << (*this)(i, j);
            }
            std::cout << std::endl;
        }
    }

    template <typename T>
    Matrix<T> Matrix<T>::copy() {
        Matrix<T> out(row, col, feat);
        return out;
    }
}


#endif //SEGM_MATRIX_H
