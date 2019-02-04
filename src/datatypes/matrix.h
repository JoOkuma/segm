
#ifndef SEGM_MATRIX_H
#define SEGM_MATRIX_H

#include "datatypes/vector.h"

#include <iostream>
#include <stdexcept>
#include <string>

#include "cblas.h"

namespace segm
{
    template<typename T>
    class Matrix
    {
    public:
        Matrix(int _row, int _col);
        Matrix(int _row, int _col, const T *array);
        Matrix(int _row, int _col, T *array, bool alloc = true);
        Matrix(const Matrix<T> &matrix);
        ~Matrix();

        Matrix<T> operator*(Matrix<T> &B) { return mult(B); };

        Matrix<T> &operator=(const Matrix<T> &matrix);
        Matrix<T> &operator+=(const Matrix<T> &matrix);
        Matrix<T> &operator-=(const Matrix<T> &matrix);

        Matrix<T> &operator=(T value);
        Matrix<T> &operator+=(T value);
        Matrix<T> &operator-=(T value);
        Matrix<T> &operator*=(T value);
        Matrix<T> &operator/=(T value);

        int getRow() const { return row; }

        int getCol() const { return col; }

        T* getFeats(int i = 0) const { return &feat[row_index[i]]; }

        void setFeats(const T *array);

        T &operator()(int _row, int _col) const { return feat[row_index[_row] + _col]; }

        T &operator()(int idx) const { return feat[idx]; }

        Matrix<T> t() const;

        void mult(const Matrix<T> &B, Matrix<T> &out, bool A_transpose = false,
                  bool B_transpose = false, float alpha = 1.0) const;

        Matrix<T> mult(const Matrix<T> &B, bool A_transpose = false,
                       bool B_transpose = false, float alpha = 1.0) const;

        void print() const;

        Matrix<T> copy() const;

        template<typename U>
        Matrix<U> convert() const;

        Matrix<T> trimCols(int begin, int end) const;
        Matrix<T> trimRows(int begin, int end) const;

        T sum() const;

    protected:
        int row;
        int col;

        T *feat = nullptr;
        int *row_index = nullptr;
        bool allocated = true;

    private:

        void checkDimensions(int _row, int _col, const char *msg);

    };

    template<typename T>
    Matrix<T>::Matrix(int _row, int _col) :
            row(_row), col(_col)
    {
        feat = new T[row * col]();
        row_index = new int[row];
        for (int i = 0, r = 0; i < row; i++, r += col)
            row_index[i] = r;

        allocated = true;
    }

    template<typename T>
    Matrix<T>::Matrix(int _row, int _col, const T *array) :
            Matrix(_row, _col) {
        setFeats(array);
    }

    template<typename T>
    Matrix<T>::Matrix(int _row, int _col, T *array, bool alloc) :
            row(_row), col(_col)
    {
        if (alloc) {
            feat = new T[row * col]();
            allocated = true;
            setFeats(array);
        } else {
            feat = array;
            allocated = false;
        }

        row_index = new int[row];
        for (int i = 0, r = 0; i < row; i++, r += col)
            row_index[i] = r;
    }

    template<typename T>
    Matrix<T>::Matrix(const Matrix<T> &matrix) :
            Matrix(matrix.getRow(), matrix.getCol()) {
        setFeats(matrix.getFeats());
    }

    template<typename T>
    Matrix<T>::~Matrix()
    {
        if (allocated)
            delete[] feat;
        delete[] row_index;
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator=(const Matrix<T> &matrix)
    {
        if (row != matrix.getRow() || col != matrix.getCol() || !allocated)
        {
            if (allocated)
                delete[] feat;
            delete[] row_index;

            row = matrix.getRow();
            col = matrix.getCol();

            feat = new T[row * col];
            row_index = new int[row];
            for (int i = 0, r = 0; i < row; i++, r += col)
                row_index[i] = r;
            allocated = true;
        }

        setFeats(matrix.getFeats());

        return (*this);
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &matrix)
    {
        checkDimensions(matrix.getRow(), matrix.getCol(),
                        "Matrix dimensions must match for addition");

        for (int i = 0; i < row * col; i++)
            feat[i] += matrix(i);

        return (*this);
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &matrix)
    {
        checkDimensions(matrix.getRow(), matrix.getCol(),
                        "Matrix dimensions must match for subtraction");

        for (int i = 0; i < row * col; i++)
            feat[i] -= matrix(i);

        return (*this);
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator=(T value) {
        for (int i = 0; i < col * row; i++) {
            feat[i] = value;
        }

        return (*this);
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator+=(T value) {
        for (int i = 0; i < col * row; i++) {
            feat[i] += value;
        }

        return (*this);
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator-=(T value) {
        for (int i = 0; i < col * row; i++) {
            feat[i] -= value;
        }

        return (*this);
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator*=(T value) {
        for (int i = 0; i < col * row; i++) {
            feat[i] *= value;
        }

        return (*this);
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator/=(T value) {
        for (int i = 0; i < col * row; i++) {
            feat[i] /= value;
        }

        return (*this);
    }

    template<typename T>
    void Matrix<T>::setFeats(const T *array) {
        for (int i = 0; i < row * col; i++)
            feat[i] = array[i];
    }

    template<typename T>
    Matrix<T> Matrix<T>::t() const {
        Matrix<T> out(col, row);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                out(j, i) = (*this)(i, j);
            }
        }
        return out;
    }

    template<typename T>
    void Matrix<T>::mult(const Matrix<T> &_B, Matrix<T> &out, bool A_transpose,
                         bool B_transpose, float alpha) const
    {
        Matrix<T> A = ((A_transpose) ? this->t() : (*this));
        Matrix<T> B = ((B_transpose) ? _B.t() : _B);

        if (A.getRow() != out.getRow() || B.getCol() != out.getCol())
            throw std::runtime_error("Dimensions of output matrix (" + std::to_string(out.getRow()) + ", " +
                                     std::to_string(out.getCol()) + ") doesn't match A (" + std::to_string(row) +
                                     ", " + std::to_string(col) + ((A_transpose) ? ")^T" : ")") + "and B (" +
                                     std::to_string(_B.getRow()) + ",  " + std::to_string(_B.getCol()) +
                                     ((B_transpose) ? ")^T" : ")"));

        if (A.getCol() != B.getRow())
            throw std::runtime_error("Cannot multiply matrices A is (" + std::to_string(row) + ", " +
                                     std::to_string(col) + ((A_transpose) ? ")^T" : ")") + "and B is (" +
                                     std::to_string(_B.getRow()) + ",  " + std::to_string(_B.getCol()) +
                                     ((B_transpose) ? ")^T" : ")"));

        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for (int i = 0; i < A.getRow(); i++) {
            for (int j = 0; j < B.getCol(); j++) {
                out(i, j) = 0;
                for (int k = 0; k < A.getCol(); k++) {
                    out(i, j) += alpha * A(i, k) * B(k, j);
                }
            }
        }
    }

    template<typename T>
    Matrix<T> Matrix<T>::mult(const Matrix<T> &B, bool A_transpose,
                                      bool B_transpose, float alpha) const
    {
        int A_row = (A_transpose) ? col : row;
        int B_col = (B_transpose) ? B.getRow() : B.getCol();

        Matrix<T> C(A_row, B_col);

        mult(B, C, A_transpose, B_transpose, alpha);

        return C;
    }

    template<typename T>
    void Matrix<T>::print() const  {
        for (int i = 0; i < row; i++) {
            std::cout << (*this)(i, 0);
            for (int j = 1; j < col; j++) {
                std::cout << " " << (*this)(i, j);
            }
            std::cout << std::endl;
        }
    }

    template<typename T>
    Matrix<T> Matrix<T>::copy() const {
        Matrix<T> out(row, col, feat);
        return out;
    }

    template<typename T>
    template<typename U>
    Matrix<U> Matrix<T>::convert() const {
        Matrix<U> out(row, col);
        for (int i = 0; i < row * col; i++)
            out(i) = (U) feat[i];
        return out;
    }


    template<typename T>
    Matrix<T> Matrix<T>::trimCols(int begin, int end) const
    {
        if (begin < 0 || end > col || (begin >= end))
            throw std::invalid_argument("Trim columns range not valid");

        Matrix<T> trimmed(row, end - begin);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < trimmed.getCol(); j++) {
                trimmed(i, j) = (*this)(i, j + begin);
            }
        }
        return trimmed;
    }

    template<typename T>
    Matrix<T> Matrix<T>::trimRows(int begin, int end) const
    {
        if (begin < 0 || end > row || (begin >= end))
            throw std::invalid_argument("Trim rows range not valid");

        Matrix<T> trimmed(end - begin, col);

        for (int i = 0; i < trimmed.getRow(); i++) {
            for (int j = 0; j < col; j++) {
                trimmed(i, j) = (*this)(i + begin, j);
            }
        }
        return trimmed;
    }

    template<typename T>
    T Matrix<T>::sum() const {
        T res = 0;
        for (int i = 0; i < row * col; i++)
            res += feat[i];
        return res;
    }

    template<typename T>
    void Matrix<T>::checkDimensions(int _row, int _col, const char *msg) {
        if (row != _row || _col != col)
            throw std::invalid_argument(msg);
    }
}


#endif //SEGM_MATRIX_H
