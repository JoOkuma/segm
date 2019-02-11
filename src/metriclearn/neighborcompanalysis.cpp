
#include "neighborcompanalysis.h"
#include "math/statistics.h"

#include <cblas.h>
#include <cmath>

using namespace segm;

NeighborCompAnalysis::NeighborCompAnalysis(MatrixXd &_data, VectorXi &_label, int output_dim,
                                           int _iterations, double _learn_rate, bool _verbose) :
        data(_data),
        label(_label),
        L(output_dim, data.cols()),
        size(static_cast<int>(data.rows())),
        d_in(static_cast<int>(data.cols())),
        d_out(output_dim)
{
    if (data.rows() != label.rows())
        throw std::invalid_argument("Number of samples on data and label must be the same");

    if (d_out > d_in)
        throw std::invalid_argument("Output dimension cannot be grater than input dimension");

    iterations = _iterations;
    learn_rate = _learn_rate;
    verbose = _verbose;
    executed = false;

    L = PCA(data, output_dim);
}


MatrixXf NeighborCompAnalysis::transform(const MatrixXf &data) const
{
    if (!executed)
        throw std::runtime_error("Neighborhood Component Analysis must"
                                 "be executed before transforming space");

    return data * L.cast<float>().transpose();
}


MatrixXd NeighborCompAnalysis::transform(const MatrixXd &data) const
{
    if (!executed)
        throw std::runtime_error("Neighborhood Component Analysis must"
                                 "be executed before transforming space");

    return data * L.transpose();
}


void NeighborCompAnalysis::train()
{

    VectorXi idx = randomIndexes(size);

    VectorXd softmax(size);

    VectorXd diff(d_in);
    VectorXd Lx_buffer(d_out);

    MatrixXd first(d_in, d_in);
    MatrixXd second(d_in, d_in);
    MatrixXd grad(d_out, d_in);

    for (int ite = 0; ite < iterations; ite++)
    {
        int i = idx[ite % size];

        double softmax_norm = 0.0;
        for (int k = 0; k < size; k++) {
            if (k == i) {
                softmax[k] = 0.0;
                continue;
            }

            diff.noalias() = data.row(i) - data.row(k);
            Lx_buffer.noalias() = L * diff;

            softmax[k] = diff.squaredNorm();
            softmax_norm += softmax[k];
        }

        softmax /= softmax_norm;

        double p_ik = 0.0;
        for (int k = 0; k < size; k++) {
            if (label[i] == label[k])
                p_ik += softmax[k];
        }

        first.setConstant(0.0);
        second.setConstant(0.0);

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

        for (int j = 0; j < d_in; j++)
            for (int k = 0; k < d_in; k++)
                first(j, k) = first(j, k) * p_ik - second(j, k);

        grad.noalias() = L * first;
        L *= 2 * learn_rate;

        L += grad;
    }
}
