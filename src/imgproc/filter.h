//
// Created by jookuma on 01/02/19.
//

#ifndef SEGM_FILTER_H
#define SEGM_FILTER_H

#include "datatypes/image.h"

namespace segm {

    class Filter {

    public:
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

        /* filters */

        /**
         * @param size          window size
         * @param sigma         standard deviation
         * @details             truncated gaussian, length of distribution tails defined by `size`, sum of weights = 1
         * @return
         */
        static Filter gaussian(int size, float sigma);

        static Image<float> convolve(Image<float> &image, Filter &filter);

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

    };
}

#endif //SEGM_FILTER_H
