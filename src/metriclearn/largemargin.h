
#ifndef SEGM_LARGEMARGIN_H
#define SEGM_LARGEMARGIN_H

#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

namespace segm
{
    class LargeMargin {

    public:
        LargeMargin(MatrixXd &data, VectorXi &label, int output_dim, int k_targets, int k_impostors = 100000,
                    int iterations = 1000, double learn_rate = 1e-7, bool verbose = false);

        template <typename T>
        Matrix<T, Dynamic, Dynamic> transform(const Matrix<T, Dynamic, Dynamic> &data) const
        { return data.lazyProduct(L.transpose()); }

        void train();

        void setLearningRate(double rate) { initial_learn_rate = rate; }

        void setIteration(int ite) { iterations = ite; }

        void setVerbose(bool _verbose) { verbose = _verbose; }

        MatrixXd getTransform() const { return L; }

    protected:

        typedef struct imp_set {
            int example;
            int impostor;
            int target;

            bool operator!=(const struct imp_set &b) const {
                return (target != b.target || impostor != b.impostor || example != b.example);
            }
            bool operator==(const struct imp_set &b) const { return !(*this != b); }

            bool operator<(const struct imp_set &b) const {
                if (example < b.example) return true;
                if (example > b.example) return false;

                if (impostor < b.impostor) return true;
                if (impostor > b.impostor) return false;

                return (target < b.target);
            }

        } impSet;

        MatrixXd &data;
        VectorXi &label;

        MatrixXd *current_data;

        MatrixXd L;

        double *dist_table = nullptr;
        bool dist_table_computed = false;

        int size;
        int d_in;
        int d_out;

        int k_targets;
        int k_impostors;

        double initial_learn_rate;
        double current_learn_rate;
        int iterations;
        bool verbose;

    private:
        void createDistTable(MatrixXd *current);
        double distance(int i, int j);

        MatrixXi findTargets();
        void checkKTargets();
        MatrixXd computeTargetGrad(const MatrixXi &target);

        virtual std::vector<impSet> findImpostors(const MatrixXi &target);
        void findDifference(const std::vector<LargeMargin::impSet> &impostors,
                            const std::vector<LargeMargin::impSet> &next_imp,
                            std::vector<LargeMargin::impSet> &missing,
                            std::vector<LargeMargin::impSet> &extra);

        MatrixXd  computeImpGrad(const std::vector<LargeMargin::impSet> &impostors);
        void computeImpGrad(const MatrixXd &imp_grad, MatrixXd &next_imp_grad,
                            const std::vector<LargeMargin::impSet> &missing,
                            const std::vector<LargeMargin::impSet> &extra);

        void updateTransform(const MatrixXd &target_grad, const MatrixXd &imp_grad,
                             MatrixXd &gradient, MatrixXd &grad_prod, MatrixXd &next_L);

        double gradLoss(const MatrixXd &grad, const MatrixXd &L_trans);


        inline double squaredNorm(const VectorXd &x, const VectorXd &y) {
            VectorXd z = x - y;
            return z.squaredNorm();
        };

    };
}


#endif //SEGM_LARGEMARGIN_H
