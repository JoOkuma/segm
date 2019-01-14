
#ifndef SEGM_LARGEMARGIN_H
#define SEGM_LARGEMARGIN_H

#include <vector>
#include "datatypes/matrix.h"

// TODO VALIDAR OVERHEAD OPEN SIMD E PARALLEL

namespace segm
{
    class LargeMargin {
    public:
        LargeMargin(int _size, int dimension, int _k_targets);

        LargeMargin(int _size, int dimension, int _k_targets, int _k_impostors);

        ~LargeMargin();

        Matrix<float> transform(const Matrix<float> &data) { return data.mult(*L, false, true); }

        void train(Matrix<float> *_data, const int *_label);

        void setLearningRate(float rate) { initial_learn_rate = rate; }

        void setIteration(int ite) { iterations = ite; }

        void setVerbose(bool _verbose) { verbose = _verbose; }

        Matrix<float> getTransform() const { return L->copy(); }

    protected:

        typedef struct imp_set {
            int example;
            int impostor;
            int target;

            bool operator!=(struct imp_set &b) {
                return (target != b.target || impostor != b.impostor || example != b.example);
            }
            bool operator==(struct imp_set &b) { return !(*this != b); }

            bool operator<(struct imp_set &b) {
                if (example < b.example) return true;
                if (example > b.example) return false;

                if (impostor < b.impostor) return true;
                if (impostor > b.impostor) return false;

                return (target < b.target);
            }

        } impSet;

        void transform(const Matrix<float> *in_data, const Matrix<float> *L_transf, Matrix<float> *out_data) {
            in_data->mult(*L_transf, *out_data, false, true);
        };

        Matrix<float> *data = nullptr;
        const int *label = nullptr;

        Matrix<float> *current_data = nullptr;

        Matrix<float> *next_Ldata = nullptr;
        Matrix<float> *L = nullptr;
        Matrix<float> *next_L = nullptr;

        float *dist_table = nullptr;
        bool dist_table_computed = false;

        Matrix<int> *target = nullptr;
        Matrix<float> *target_grad = nullptr;

        std::vector<impSet> *impostors = nullptr;
        std::vector<impSet> *next_imp = nullptr;
        std::vector<impSet> *missing = nullptr;
        std::vector<impSet> *extra = nullptr;

        Matrix<float> *imp_grad = nullptr;
        Matrix<float> *next_imp_grad = nullptr;

        Matrix<float> *gradient = nullptr;

        int size;
        int d;

        int k_targets;
        int k_impostors = 100000;

        float initial_learn_rate = 1e-7;
        float current_learn_rate;
        int iterations = 1000;
        bool verbose = false;

    private:
        void allocAuxiliary();
        void freeAuxiliary();

        void createDistTable(Matrix<float> *current);
        float distance(int i, int j);

        void findTargets();
        void checkKTargets();
        void computeTargetGrad();

        virtual std::vector<impSet> *findImpostors();
        void findDifference();
        void computeImpGrad();

        void computeJointGrad();
        void updateTransform();

        float gradLoss(const Matrix<float> *grad, const Matrix<float> *L_trans);

        inline float squaredl2norm(const float *x, const float *y)
        {
            float dist = 0.0;
            for (int i = 0; i < d; i++) {
                float diff = x[i] - y[i];
                dist += diff * diff;
            }
            return dist;
        };

        void swap(Matrix<float> **A, Matrix<float> **B) {
            Matrix<float> *tmp = *A;
            *A = *B;
            *B = tmp;
        }
    };
}


#endif //SEGM_LARGEMARGIN_H
