
#include "neighborcompanalysis.h"
#include "math/statistics.h"

#include <cblas.h>
#include <cmath>

using namespace segm;

NeighborCompAnalysis::NeighborCompAnalysis(Matrix<float> &_data, Vector<int> &_label, int output_dim,
                                           int _iterations, double _learn_rate, bool _verbose) :
        data(_data),
        label(_label),
        L(output_dim, data.getCol()),
        size(data.getRow()),
        d_in(data.getCol()),
        d_out(output_dim)
{
    if (data.getRow() != label.getSize())
        throw std::invalid_argument("Number of samples on data and label must be the same");

    if (d_out > d_in)
        throw std::invalid_argument("Output dimension cannot be grater than input dimension");

    iterations = _iterations;
    learn_rate = _learn_rate;
    verbose = _verbose;
    executed = false;

    Matrix<double> tmp_data = data.convert<double>();
    L = PCA(tmp_data, output_dim);
}


Matrix<float> NeighborCompAnalysis::transform(const Matrix<float> &data) const
{
    if (!executed)
        throw std::runtime_error("Neighborhood Component Analysis must"
                                 "be executed before transforming space");

    return data.mult(L.convert<float>(), false, true);
}


Matrix<double> NeighborCompAnalysis::transform(const Matrix<double> &data) const
{
    if (!executed)
        throw std::runtime_error("Neighborhood Component Analysis must"
                                 "be executed before transforming space");

    return data.mult(L, false, true);
}


void NeighborCompAnalysis::train()
{

    Vector<int> idx = randomIndexes(size);

    Vector<double> softmax(size);

    Vector<double> diff(d_in);
    Vector<double> Lx_buffer(d_out);

    Matrix<double> first(d_in, d_in);
    Matrix<double> second(d_in, d_in);
    Matrix<double> grad(d_out, d_in);

    for (int ite = 0; ite < iterations; ite++)
    {
        int i = idx[ite % size];

        double softmax_norm = 0.0;
        for (int k = 0; k < size; k++) {
            if (k == i) {
                softmax[k] = 0.0;
                continue;
            }

            for (int j = 0; j < d_in; j++) {
                diff[j] = data(i, j) - data(k, j);
            }

            softmax[k] = transSqrL2Norm(diff, L, Lx_buffer);
            softmax_norm += softmax[k];
        }

        softmax /= softmax_norm;

        double p_ik = 0.0;
        for (int k = 0; k < size; k++) {
            if (label[i] == label[k])
                p_ik += softmax[k];
        }

        first = 0.0;
        second = 0.0;

        for (int di = 0; di < d_in; di++) {
            for (int dj = 0; dj < d_in; dj++) {
                for (int k = 0; k < size; k++)
                {
                    double x = + data(i, di) * data(i, dj)
                               + data(k, di) * data(k, dj)
                               - data(i, di) * data(k, dj)
                               - data(k, di) * data(i, dj);

                    x *= softmax[k];
                    first(i, dj) += x;
                    if (label[k] == label[i]) {
                        second(i, dj) += x;
                    }
                }
            }
        }

        for (int j = 0; j < d_in * d_in; j++)
            first(j) = first(j) * p_ik - second(j);

        L.mult(first, grad, false, false, 2 * learn_rate);

        L += grad;
    }
}


double NeighborCompAnalysis::transSqrL2Norm(Vector<double> &array, Matrix<double> &L, Vector<double> &buffer)
{
    double norm = 0.0;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, d_out, d_in, 1.0, array.getFeats(),
                d_in, L.getFeats(), d_in, 0.0, buffer.getFeats(), d_out);

    for (int i = 0; i < d_out; i++) {
        norm += buffer[i] * buffer[i];
    }

    return norm;
}

