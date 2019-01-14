
#include <cmath>
#include <iostream>
#include <limits>

#include "largemargin.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace segm;

LargeMargin::LargeMargin(int _size, int dimension, int _k_targets) :
    size(_size), d(dimension), k_targets(_k_targets)
{
    L = new Matrix<float>(d, d);

    for (int i = 0; i < d; i++) {
        (*L)(i, i) = 1.0;
    }
}

LargeMargin::LargeMargin(int _size, int dimension, int _k_targets, int _k_impostors) :
        size(_size), d(dimension), k_targets(_k_targets), k_impostors(_k_impostors)
{
    L = new Matrix<float>(d, d);

    for (int i = 0; i < d; i++) {
        (*L)(i, i) = 1.0;
    }
}

LargeMargin::~LargeMargin()
{
    delete L;
}

void LargeMargin::train(Matrix<float> *_data, const int *_label)
{
    data = _data;
    label = _label;
    allocAuxiliary();

    checkKTargets();
    createDistTable(data);

    if (verbose && dist_table_computed)
        std::cout << "Distance matrix computed" << std::endl;

    findTargets();
    computeTargetGrad();

    if (verbose)
        std::cout << "Target neighbours found and target static gradient computed" << std::endl;

    for (int i = 0; i < d * d; i++) {
        (*next_L)(i) = (*L)(i);
    }

    const float epsilon = 1e-7;
    const float min_learn_rate = 1e-22;

    current_learn_rate = initial_learn_rate;
    int ite = 0;
    float delta;
    float loss = std::numeric_limits<float>::max();
    do {
        ite++;
        do {
            updateTransform();
            transform(data, next_L, next_Ldata);
            createDistTable(next_Ldata);
            next_imp = findImpostors();
            findDifference();
            computeImpGrad();
            float next_loss = gradLoss(target_grad, next_L) + gradLoss(next_imp_grad, next_L) + next_imp->size();
            delta = next_loss - loss;
            if (delta > 0.0) {
                delete next_imp;
                current_learn_rate /= 2;
            } else {
                current_learn_rate *= 1.01;
                loss = next_loss;
                swap(&L, &next_L);
                swap(&imp_grad, &next_imp_grad);
                delete impostors;
                impostors = next_imp;
            }
        } while (current_learn_rate > min_learn_rate && delta > 0.0);

        if (verbose)
            std::cout << "Iteration: " << ite << std::endl
                      << "Loss: " << loss << std::endl
                      << "Delta: " << delta << std::endl
                      << "Active Impostors: " << impostors->size() << std::endl
                      << "Learning Rate: " << current_learn_rate << std::endl;

    } while (ite < iterations && !impostors->empty()
             && fabsf(delta) > epsilon && current_learn_rate > min_learn_rate);

    delete[] dist_table;
    freeAuxiliary();
}

void LargeMargin::allocAuxiliary()
{
    next_L = new Matrix<float>(d, d);
    next_Ldata = new Matrix<float>(size, d);

    gradient = new Matrix<float>(d, d);
    grad_prod = new Matrix<float>(d, d);

    target = new Matrix<int>(size, k_targets);
    target_grad = new Matrix<float>(d, d);

    impostors = new std::vector<LargeMargin::impSet>;
    missing = new std::vector<LargeMargin::impSet>;
    extra = new std::vector<LargeMargin::impSet>;

    imp_grad = new Matrix<float>(d, d);
    next_imp_grad = new Matrix<float>(d, d);
}

void LargeMargin::freeAuxiliary()
{
    delete next_L;
    delete next_Ldata;

    delete gradient;
    delete grad_prod;

    delete target;
    delete target_grad;

    delete impostors;
    delete missing;
    delete extra;

    delete imp_grad;
    delete next_imp_grad;
}

void LargeMargin::createDistTable(Matrix<float> *current)
{
    current_data = current;

    if (dist_table_computed) {
        delete[] dist_table;
        dist_table_computed = false;
    }

    dist_table = new (std::nothrow) float[(size * (size - 1)) / 2];
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
            dist_table[row_id + j] = squaredl2norm(current->getFeat(i), current->getFeat(j));
        }
    }

    if (verbose)
        std::cout << "Computed triangular distance matrix" << std::endl;

    dist_table_computed = true;
}

float LargeMargin::distance(int i, int j)
{
    // very unlikely to happen
    if (i == j) return std::numeric_limits<float>::max();

    if (!dist_table_computed)
        return squaredl2norm(current_data->getFeat(i), current_data->getFeat(j));

    int greater = (i > j) ? i : j;
    int lower = (i > j) ? j : i;
    int row_id = (greater * (greater - 1) / 2);
    return dist_table[row_id + lower];
}


void LargeMargin::findTargets()
{
    auto *max_dist = new float[size];
    auto *target_dist = new Matrix<float>(size, k_targets);

    for (int i = 0; i < k_targets * size; i++) {
        (*target_dist)(i) = std::numeric_limits<float>::max();
        (*target)(i) = -1;
    }

    for (int i = 0; i < size; i++) {
        max_dist[i] = std::numeric_limits<float>::max();
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (label[i] == label[j] && i != j)
            {
                float dist = distance(i, j);
                if (dist < max_dist[i])
                {
                    int last = k_targets - 1;
                    int idx = last;
                    while (idx > 0 && dist < (*target_dist)(i, idx - 1)) {
                        (*target_dist)(i, idx) = (*target_dist)(i, idx - 1);
                        (*target)(i, idx) = (*target)(i, idx - 1);
                        idx--;
                    }
                    (*target_dist)(i, idx) = dist;
                    (*target)(i, idx) = j;
                    max_dist[i] = (*target_dist)(i, last);
                }
            }
        }
    }

    delete[] max_dist;
    delete target_dist;
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

void LargeMargin::computeTargetGrad()
{
    for (int di = 0; di < d; di++) {
        for (int dj = 0; dj < d; dj++) {
            for (int i = 0; i < size; i++) {
                for (int k = 0; k < k_targets; k++)
                {
                    (*target_grad)(di, dj) += + (*data)(i, di) * (*data)(i, dj)
                                              + (*data)(k, di) * (*data)(k, dj)
                                              - (*data)(i, di) * (*data)(k, dj)
                                              - (*data)(k, di) * (*data)(i, dj);
                }
            }
        }
    }
}

std::vector<LargeMargin::impSet> *LargeMargin::findImpostors()
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
                float imp_dist = distance(i, ik);
                bool is_imp = false;
                for (int j = 0; j < k_targets; j++)
                {
                    int tk = (*target)(i, j);
                    float diff = distance(i, tk) + 1.0f - imp_dist;
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

    auto *out = new std::vector<LargeMargin::impSet>;
    for (unsigned long i = 0; i < n_threads; i++) {
        out->insert(out->end(), imp_arr[i].begin(), imp_arr[i].end());
    }

    return out;
}

void LargeMargin::computeImpGrad()
{
    for (int i = 0; i < d * d; i++)
        (*next_imp_grad)(i) = (*imp_grad)(i);

    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (int di = 0; di < d; di++) {
        for (int dj = 0; dj < d; dj++)
        {
            /* adding */
            for (impSet s : *missing) {
                int e = s.example;
                int ik = s.impostor;
                int tk = s.target;

                /* impostors */
                (*next_imp_grad)(di, dj) += - (*data)(e, di) * (*data)(e, dj)
                                            - (*data)(ik, di) * (*data)(ik, dj)
                                            + (*data)(e, di) * (*data)(ik, dj)
                                            + (*data)(ik, di) * (*data)(e, dj);

                /* target */
                (*next_imp_grad)(di, dj) += + (*data)(e, di) * (*data)(e, dj)
                                            + (*data)(tk, di) * (*data)(tk, dj)
                                            - (*data)(e, di) * (*data)(tk, dj)
                                            - (*data)(tk, di) * (*data)(e, dj);
            }

            /* subtracting */
            for (impSet s : *extra) {
                int e = s.example;
                int ik = s.impostor;
                int tk = s.target;

                /* impostors */
                (*next_imp_grad)(di, dj) -= - (*data)(e, di) * (*data)(e, dj)
                                            - (*data)(ik, di) * (*data)(ik, dj)
                                            + (*data)(e, di) * (*data)(ik, dj)
                                            + (*data)(ik, di) * (*data)(e, dj);

                /* target */
                (*next_imp_grad)(di, dj) -= + (*data)(e, di) * (*data)(e, dj)
                                            + (*data)(tk, di) * (*data)(tk, dj)
                                            - (*data)(e, di) * (*data)(tk, dj)
                                            - (*data)(tk, di) * (*data)(e, dj);

            }
        }
    }
}

/* find extra and missing impostors from `next_imp to `impostors` */
void LargeMargin::findDifference()
{
    unsigned long smaller_size = (next_imp->size() < impostors->size()) ?
                                 next_imp->size() : impostors->size();
    missing->clear();
    extra->clear();
    unsigned long i = 0, j = 0;
    for (i = 0, j = 0; i < smaller_size && j < smaller_size; i++, j++) {
        if (next_imp->at(i) != impostors->at(j)) {
            if (next_imp->at(i) < impostors->at(j)) {
                missing->push_back(next_imp->at(i));
                j--;
            } else {
                extra->push_back(impostors->at(j));
                i--;
            }
        }
    }

    for (; i < next_imp->size(); i++) {
        missing->push_back(next_imp->at(i));
    }
    for (; j < impostors->size(); j++) {
        extra->push_back(impostors->at(j));
    }
}

void LargeMargin::computeJointGrad()
{
    for (int i = 0; i < d * d; i++) {
        (*gradient)(i) = (*target_grad)(i) + (*imp_grad)(i);
    }
}

void LargeMargin::updateTransform()
{
    computeJointGrad();

    L->mult(*gradient, *grad_prod);
    for (int i = 0; i < d * d; i++) {
        (*next_L)(i) = (*L)(i) - 2 * current_learn_rate * (*grad_prod)(i);
    }
}

float LargeMargin::gradLoss(const Matrix<float> *grad, const Matrix<float> *L_trans)
{
    Matrix<float> M = L_trans->mult(*L_trans, true, false);

    float loss = 0.0f;
    for (int i = 0; i < d * d; i++) {
        loss += M(i) * (*grad)(i);
    }

    return loss;
}

