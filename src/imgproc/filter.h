//
// Created by jookuma on 01/02/19.
//

#ifndef SEGM_FILTER_H
#define SEGM_FILTER_H

#include "datatypes/image.h"

#include <cmath>

namespace segm {

    class Filter {

    public:

        enum ConductionFunction {
            exponential,
            quadratic,
        };

        Filter(int width, int height) : Filter(width, height, 1) { }
        Filter(int width, int height, int band);
        ~Filter();

        int getWidth() const { return w; }
        int getHeight() const { return h; }
        int getBands() const { return b; }

        float &operator()(int x, int y, int b) { return feat[row[y] + col[x] + b]; }
        /* only for 1 band filters */
        float &operator()(int x, int y) { return feat[row[y] + x]; }
        float &operator()(int p) { return feat[p]; }

        static Image<float> convolve(Image<float> &image, Filter &filter);

        /* filters */

        /**
         * @param size          window size
         * @param sigma         standard deviation
         * @details             truncated gaussian, length of distribution tails defined by `size`, sum of weights = 1
         * @return
         */
        static Filter gaussian(int size, float sigma);

        static Image<float> anisotropic(const Image<float> &image, ConductionFunction fun_type,
                                        int iters = 50, float lambda = 0.1f, float kappa = 0.5f);

    protected:
        int h = 0; /* h */
        int w = 0; /* w */
        int b = 0; /* bands */
        float *feat = nullptr;

        int *row = nullptr; /* w plus band length padding */
        int *col = nullptr; /* band padding only */
        int *row_index = nullptr; /* w padding only */

        int *shift_h = nullptr;
        int *shift_w = nullptr;

    private:

        static inline double expConduc(const float *f1, const float *f2, int b, float kappa)
        {
            double dot_product = 0.0;
            for (int i = 0; i < b; i++) {
                double diff = f1[b] - f2[b];
                dot_product += diff * diff;
            }
            return 1.0 / (1.0 + dot_product / (kappa * kappa));
        }


        static inline double quadConduct(const float *f1, const float *f2, int b, float kappa)
        {
            double dot_product = 0.0;
            for (int i = 0; i < b; i++) {
                double diff = f1[b] - f2[b];
                dot_product += diff * diff;
            }
            return exp(-dot_product / (kappa * kappa));
        }


    };
}

#endif //SEGM_FILTER_H
