
#include <cmath>
#include <iostream>
#include <limits>

#include "metriclearn/largemargin.h"
#include "math/statistics.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace segm;

LargeMargin::LargeMargin(const Matrix<float> &_data, Vector<int> &_label, int output_dim, int _k_targets, int _k_impostors,
                         int _iterations, double _learn_rate, bool _verbose) :
        label(_label), L(output_dim, _data.getCol()), size(_data.getRow()), d_in(_data.getCol()), d_out(output_dim),
        k_targets(_k_targets), k_impostors(_k_impostors), initial_learn_rate(_learn_rate),
        iterations(_iterations), verbose(_verbose)
{

    checkKTargets();
}


LargeMargin::~LargeMargin()
{
}

void LargeMargin::train(const Matrix<float> &_data)
{
    data = new Matrix<double>(_data.convert<double>());
    L = Matrix<double>(PCA(*data, d_out));

    createDistTable(data);

    if (verbose && dist_table_computed)
        std::cout << "Distance matrix computed" << std::endl;

    Matrix<int> target = findTargets();
    Matrix<double> target_grad = computeTargetGrad(target);

    if (verbose)
        std::cout << "Target neighbours found and target static gradient computed" << std::endl;

    Matrix<double> next_L(L);
    Matrix<double> next_Ldata = data->mult(L, false, true);

    createDistTable(&next_Ldata);

    std::vector<LargeMargin::impSet> impostors = findImpostors(target);
    Matrix<double> imp_grad = computeImpGrad(impostors);

    /* buffers */
    Matrix<double> next_imp_grad(d_in, d_in);
    Matrix<double> gradient(d_in, d_in);
    Matrix<double> grad_prod(d_out, d_in);

    const double epsilon = 1e-7;
    const double min_learn_rate = 1e-22;

    current_learn_rate = initial_learn_rate;
    int ite = 0;
    double delta;
    double loss = gradLoss(target_grad, next_L) + gradLoss(imp_grad, next_L) + impostors.size();
    do {
        ite++;
        do {
            updateTransform(target_grad, imp_grad, gradient, grad_prod, next_L);
            transform(*data, next_L, next_Ldata);
            createDistTable(&next_Ldata);

            std::vector<LargeMargin::impSet> next_imp = findImpostors(target);
            std::vector<LargeMargin::impSet> missing;
            std::vector<LargeMargin::impSet> extra;
            findDifference(impostors, next_imp, missing, extra);

            computeImpGrad(imp_grad, next_imp_grad, missing, extra);

            double next_loss = gradLoss(target_grad, next_L) + gradLoss(next_imp_grad, next_L) + next_imp.size();

            delta = next_loss - loss;
            if (delta > 0.0) {
                current_learn_rate /= 2;
            } else {
                current_learn_rate *= 1.01;
                loss = next_loss;
                L.swap(next_L);
                imp_grad.swap(next_imp_grad);
                impostors = next_imp;
            }
        } while (current_learn_rate > min_learn_rate && delta > 0.0);

        if (verbose)
            std::cout << "Iteration: " << ite << std::endl
                      << "Loss: " << loss << std::endl
                      << "Delta: " << delta << std::endl
                      << "Active Impostors: " << impostors.size() << std::endl
                      << "Learning Rate: " << current_learn_rate << std::endl;

    } while (ite < iterations && !impostors.empty()
             && fabs(delta) > epsilon && current_learn_rate > min_learn_rate);

    delete[] dist_table;
}


void LargeMargin::createDistTable(Matrix<double> *current)
{
    current_data = current;

    if (dist_table_computed) {
        delete[] dist_table;
        dist_table_computed = false;
    }

    dist_table = new (std::nothrow) double[(size * (size - 1)) / 2];
    if (!dist_table) {
        dist_table_computed = false;
        return;
    }

    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (int i = 1; i < size; i++) {
        int row_id = (i * (i - 1)) / 2;
        for (int j = 0; j < i; j++) {
            dist_table[row_id + j] = squaredl2norm(current->getFeats(i), current->getFeats(j), current->getCol());
        }
    }

    if (verbose)
        std::cout << "Computed triangular distance matrix" << std::endl;

    dist_table_computed = true;
}

double LargeMargin::distance(int i, int j)
{
    // very unlikely to happen
    if (i == j) return std::numeric_limits<double>::max();

    if (!dist_table_computed)
        return squaredl2norm(current_data->getFeats(i), current_data->getFeats(j), current_data->getCol());

    int greater = (i > j) ? i : j;
    int lower = (i > j) ? j : i;
    int row_id = (greater * (greater - 1) / 2);
    return dist_table[row_id + lower];
}


Matrix<int> LargeMargin::findTargets()
{
    Vector<double> max_dist(size);
    Matrix<int> target(size, k_targets);
    Matrix<double> target_dist(size, k_targets);

    target_dist = std::numeric_limits<double>::max();
    target = -1;
    max_dist = std::numeric_limits<double>::max();

    int last = k_targets - 1;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (label[i] == label[j] && i != j)
            {
                double dist = distance(i, j);
                if (dist < max_dist[i])
                {
                    int idx = last;
                    while (idx > 0 && dist < target_dist(i, idx - 1)) {
                        target_dist(i, idx) = target_dist(i, idx - 1);
                        target(i, idx) = target(i, idx - 1);
                        idx--;
                    }
                    target_dist(i, idx) = dist;
                    target(i, idx) = j;
                    max_dist[i] = target_dist(i, last);
                }
            }
        }
    }

    return target;
}

void LargeMargin::checkKTargets()
{
    int max_lab = -1;
    for (int i = 0; i < size; i++) {
        if (max_lab < label[i])
            max_lab = label[i];
    }
    max_lab += 1;

    int *lab_count = new int[max_lab]();
    for (int i = 0; i < size; i++) {
        lab_count[label[i]] += 1;
    }

    int min = std::numeric_limits<int>::max();
    for (int i = 0; i < max_lab; i++) {
        if (min > lab_count[i])
            min = lab_count[i];
    }
    min -=1;

    if (min <  k_targets) {
        k_targets = min;
        std::cout << "Number of classes reduced to " << min << ", number of samples"
                "with same class not sufficient." << std::endl;
    }

    delete[] lab_count;
}

Matrix<double> LargeMargin::computeTargetGrad(const Matrix<int> &target)
{
    Matrix<double> target_grad(d_in, d_in);

    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (int di = 0; di < d_in; di++) {
        for (int dj = 0; dj < d_in; dj++) {
            for (int i = 0; i < size; i++) {
                for (int k = 0; k < k_targets; k++)
                {
                    int tk = target(i, k);
                    target_grad(di, dj) += + (*data)(i, di) *  (*data)(i, dj)
                                           + (*data)(tk, di) * (*data)(tk, dj)
                                           - (*data)(i, di) *  (*data)(tk, dj)
                                           - (*data)(tk, di) * (*data)(i, dj);
                }
            }
        }
    }

    return target_grad;
}

std::vector<LargeMargin::impSet> LargeMargin::findImpostors(const Matrix<int> &target)
{
    #ifndef _OPENMP
    unsigned long n_threads = 1;
    std::vector<std::vector<LargeMargin::impSet>> imp_arr(1);
    #else
    unsigned long  n_threads = 8;
    std::vector<std::vector<LargeMargin::impSet>> imp_arr(n_threads);
    #pragma omp parallel for num_threads(n_threads)
    #endif
    for (int i = 0; i < size; i++)
    {
        #ifdef _OPENMP
        int thread_id = omp_get_thread_num();
        #else
        int thread_id = 0;
        #endif
        int count = 0;
        for (int ik = 0; ik < size; ik++)
        {
            if (count >= k_impostors) break;
            if (label[i] != label[ik])
            {
                double imp_dist = distance(i, ik);
                bool is_imp = false;
                for (int j = 0; j < k_targets; j++)
                {
                    int tk = target(i, j);
                    double diff = distance(i, tk) + 1.0 - imp_dist;
                    if (diff > 0) {
                        LargeMargin::impSet s = {.example = i, .impostor = ik, .target = tk};
                        imp_arr[thread_id].push_back(s);
                        is_imp = true;
                    }
                }
                if (is_imp) count++;
            }
        }
    }

    std::vector<LargeMargin::impSet> out;
    for (unsigned long i = 0; i < n_threads; i++) {
        out.insert(out.end(), imp_arr[i].begin(), imp_arr[i].end());
    }

    return out;
}

/* find extra and missing impostors from `next_imp to `impostors` */
void LargeMargin::findDifference(const std::vector<LargeMargin::impSet> &impostors,
                                 const std::vector<LargeMargin::impSet> &next_imp,
                                 std::vector<LargeMargin::impSet> &missing,
                                 std::vector<LargeMargin::impSet> &extra)
{
    unsigned long smaller_size = (next_imp.size() < impostors.size()) ?
                                 next_imp.size() : impostors.size();

    missing.clear();
    extra.clear();
    unsigned long i = 0, j = 0;
    for (i = 0, j = 0; i < smaller_size && j < smaller_size; i++, j++) {
        if (next_imp.at(i) != impostors.at(j)) {
            if (next_imp.at(i) < impostors.at(j)) {
                missing.push_back(next_imp.at(i));
                j--;
            } else {
                extra.push_back(impostors.at(j));
                i--;
            }
        }
    }

    for (; i < next_imp.size(); i++) {
        missing.push_back(next_imp.at(i));
    }
    for (; j < impostors.size(); j++) {
        extra.push_back(impostors.at(j));
    }
}

Matrix<double> LargeMargin::computeImpGrad(const std::vector<LargeMargin::impSet> &impostors)
{
    Matrix<double> imp_grad(d_in, d_in);

    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (int di = 0; di < d_in; di++) {
        for (int dj = 0; dj < d_in; dj++)
        {
            /* adding */
            for (impSet s : impostors) {
                int e = s.example;
                int ik = s.impostor;
                int tk = s.target;

                /* impostors */
                imp_grad(di, dj) += - (*data)(e, di) * (*data)(e, dj)
                                    - (*data)(ik, di) * (*data)(ik, dj)
                                    + (*data)(e, di) * (*data)(ik, dj)
                                    + (*data)(ik, di) * (*data)(e, dj);

                /* target */
                imp_grad(di, dj) += + (*data)(e, di) * (*data)(e, dj)
                                    + (*data)(tk, di) * (*data)(tk, dj)
                                    - (*data)(e, di) * (*data)(tk, dj)
                                    - (*data)(tk, di) * (*data)(e, dj);
            }
        }
    }

    return imp_grad;
}


void LargeMargin::computeImpGrad(const Matrix<double> &imp_grad, Matrix<double> &next_imp_grad,
                                 const std::vector<LargeMargin::impSet> &missing,
                                 const std::vector<LargeMargin::impSet> &extra)
{
    next_imp_grad = imp_grad;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int di = 0; di < d_in; di++) {
        for (int dj = 0; dj < d_in; dj++)
        {
            /* adding */
            for (impSet s : missing) {
                int e = s.example;
                int ik = s.impostor;
                int tk = s.target;

                /* impostors */
                next_imp_grad(di, dj) += - (*data)(e, di) * (*data)(e, dj)
                                         - (*data)(ik, di) * (*data)(ik, dj)
                                         + (*data)(e, di) * (*data)(ik, dj)
                                         + (*data)(ik, di) * (*data)(e, dj);

                /* target */
                next_imp_grad(di, dj) += + (*data)(e, di) * (*data)(e, dj)
                                         + (*data)(tk, di) * (*data)(tk, dj)
                                         - (*data)(e, di) * (*data)(tk, dj)
                                         - (*data)(tk, di) * (*data)(e, dj);
            }

            /* subtracting */
            for (impSet s : extra) {
                int e = s.example;
                int ik = s.impostor;
                int tk = s.target;

                /* impostors */
                next_imp_grad(di, dj) -= - (*data)(e, di) * (*data)(e, dj)
                                         - (*data)(ik, di) * (*data)(ik, dj)
                                         + (*data)(e, di) * (*data)(ik, dj)
                                         + (*data)(ik, di) * (*data)(e, dj);

                /* target */
                next_imp_grad(di, dj) -= + (*data)(e, di) * (*data)(e, dj)
                                         + (*data)(tk, di) * (*data)(tk, dj)
                                         - (*data)(e, di) * (*data)(tk, dj)
                                         - (*data)(tk, di) * (*data)(e, dj);

            }
        }
    }
}


void LargeMargin::updateTransform(const Matrix<double> &target_grad, const Matrix<double> &imp_grad,
                                  Matrix<double> &gradient, Matrix<double> &grad_prod, Matrix<double> &next_L)
{
    for (int i = 0; i < d_in * d_in; i++) {
        gradient(i) = target_grad(i) + imp_grad(i);
    }

    L.mult(gradient, grad_prod);
    for (int i = 0; i < d_in * d_out; i++) {
        next_L(i) = L(i) - 2 * current_learn_rate * grad_prod(i);
    }
}

double LargeMargin::gradLoss(const Matrix<double> &grad, const Matrix<double> &L_trans)
{
    Matrix<double> M = L_trans.mult(L_trans, true, false);

    double loss = 0.0f;
    for (int i = 0; i < d_in * d_in; i++) {
        loss += M(i) * grad(i);
    }

    return loss;
}

