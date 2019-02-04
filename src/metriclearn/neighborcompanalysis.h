
#ifndef SEGM_NEIGHCOMPANALYSIS_H
#define SEGM_NEIGHCOMPANALYSIS_H

#include "datatypes/matrix.h"

namespace segm {

    class NeighborCompAnalysis {

    public:

        NeighborCompAnalysis(Matrix<float> &data, Vector<int> &label, int output_dim, int _iterations = 10000,
                             float _learn_rate = 1e-2, bool _verbose = false);

        Matrix<float> transform(const Matrix<float> &data) const;
        Matrix<double> transform(const Matrix<double> &data) const;

        void train();

    private:

        Matrix<float> &data;
        Vector<int> &label;
        Matrix<double> L;

        int size;       // input data length
        int d_in;       // input dimension
        int d_out;      // output dimension

        float learn_rate = 1e-2;
        int iterations = 10000;
        bool verbose = false;
        bool executed = false;

        double transSqrL2Norm(Vector<double> &array, Matrix<double> &L, Vector<double> &buffer);

    };
}


#endif //SEGM_NEIGHCOMPANALYSIS_H
