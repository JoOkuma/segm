
#ifndef SEGM_LARGEMARGIN_H
#define SEGM_LARGEMARGIN_H

#include <vector>
#include "datatypes/matrix.h"

namespace segm
{
    class LargeMargin {

    public:
        LargeMargin(const Matrix<float> &data, Vector<int> &label, int output_dim, int k_targets, int k_impostors = 100000,
                    int iterations = 1000, double learn_rate = 1e-7, bool verbose = false);

        ~LargeMargin();

        Matrix<double> transform(const Matrix<double> &data) const { return data.mult(L, false, true); }
        Matrix<float> transform(const Matrix<float> &data) const { return data.mult(L.convert<float>(), false, true); }

        void train(const Matrix<float> &data);

        void setLearningRate(float rate) { initial_learn_rate = rate; }

        void setIteration(int ite) { iterations = ite; }

        void setVerbose(bool _verbose) { verbose = _verbose; }

        Matrix<double> getTransform() const { return L.copy(); }

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

        void transform(const Matrix<double> &in_data,
                       const Matrix<double> &L_transf,
                       Matrix<double> &out_data) {
            in_data.mult(L_transf, out_data, false, true);
        };

        Matrix<double> *data;
        Vector<int> &label;

        Matrix<double> *current_data;

        Matrix<double> L;

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
        void createDistTable(Matrix<double> *current);
        double distance(int i, int j);

        Matrix<int> findTargets();
        void checkKTargets();
        Matrix<double> computeTargetGrad(const Matrix<int> &target);

        virtual std::vector<impSet> findImpostors(const Matrix<int> &target);
        void findDifference(const std::vector<LargeMargin::impSet> &impostors,
                            const std::vector<LargeMargin::impSet> &next_imp,
                            std::vector<LargeMargin::impSet> &missing,
                            std::vector<LargeMargin::impSet> &extra);

        Matrix<double>  computeImpGrad(const std::vector<LargeMargin::impSet> &impostors);
        void computeImpGrad(const Matrix<double> &imp_grad, Matrix<double> &next_imp_grad,
                            const std::vector<LargeMargin::impSet> &missing,
                            const std::vector<LargeMargin::impSet> &extra);

        void updateTransform(const Matrix<double> &target_grad, const Matrix<double> &imp_grad,
                             Matrix<double> &gradient, Matrix<double> &grad_prod, Matrix<double> &next_L);

        double gradLoss(const Matrix<double> &grad, const Matrix<double> &L_trans);

        inline double squaredl2norm(const double *x, const double *y, int d)
        {
            double dist = 0.0;
            for (int i = 0; i < d; i++) {
                double diff = x[i] - y[i];
                dist += diff * diff;
            }
            return dist;
        };

    };
}


#endif //SEGM_LARGEMARGIN_H
