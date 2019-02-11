
#ifndef SEGM_NEIGHCOMPANALYSIS_H
#define SEGM_NEIGHCOMPANALYSIS_H

#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace segm {

    class NeighborCompAnalysis {

    public:

        NeighborCompAnalysis(MatrixXd &data, VectorXi &label, int output_dim, int _iterations = 10000,
                             double _learn_rate = 1e-2, bool _verbose = false);

        MatrixXf transform(const MatrixXf &data) const;
        MatrixXd transform(const MatrixXd &data) const;

        void train();

    private:

        MatrixXd &data;
        VectorXi &label;
        MatrixXd L;

        int size;       // input data length
        int d_in;       // input dimension
        int d_out;      // output dimension

        double learn_rate = 1e-2;
        int iterations = 10000;
        bool verbose = false;
        bool executed = false;

    };
}


#endif //SEGM_NEIGHCOMPANALYSIS_H
