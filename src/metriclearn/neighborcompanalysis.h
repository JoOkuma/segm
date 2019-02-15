
#ifndef SEGM_NEIGHCOMPANALYSIS_H
#define SEGM_NEIGHCOMPANALYSIS_H

#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace segm {

    class NeighborCompAnalysis {

    public:

        NeighborCompAnalysis(MatrixXd &data, VectorXi &label, int output_dim, int _iterations = 10000,
                             double _learn_rate = 1e-2);

        MatrixXf transform(const MatrixXf &data) const;
        MatrixXd transform(const MatrixXd &data) const;

        void train();

        void setLearningRate(double rate) { learn_rate = rate; }
        void setIteration(int ite) { iterations = ite; }

        MatrixXd getTransform() const { return L; }

    private:

        MatrixXd &data;
        VectorXi &label;
        MatrixXd L;

        int size;       // input data length
        int d_in;       // input dimension
        int d_out;      // output dimension

        double learn_rate = 1e-2;
        int iterations = 10000;
        bool executed = false;

    };
}


#endif //SEGM_NEIGHCOMPANALYSIS_H
