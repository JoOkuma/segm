#include "filter.h"


using namespace segm;


Filter::Filter(int width, int height, int band)
{
    h = height;
    w = width;
    b = band;
    feat = new float[h * w * b];

    row = new int[h];
    row_index = new int[h];
    col = new int[w];

    shift_w = new int[w];
    shift_h = new int[h];

    int begin_w = - w / 2;
    for (int i = 0, c = 0; i < w; i++, c += b) {
        shift_w[i] = i + begin_w;
        col[i] = c;
    }

    int begin_h = - h / 2;
    for (int i = 0, r = 0; i < h; i++, r += w) {
        shift_h[i] = i + begin_h;
        row[i] = r * b;
        row_index[i] = r;
    }
}


Filter::~Filter()
{
    delete[] feat;
    delete[] row;
    delete[] row_index;
    delete[] col;
    delete[] shift_w;
    delete[] shift_h;
}


Image<float> Filter::convolve(Image<float> &image, Filter &filter)
{
    int width = image.getWidth();
    int height = image.getHeight();
    int bands = image.getBands();

    Image<float> out(width, height, bands);

    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int b = 0; b < bands; b++) {
                int b_row = (b % filter.b) * filter.h * filter.w;
                int f = 0;
                for (int x = 0; x < filter.w; x++) {
                    for (int y = 0; y < filter.h; y++) {
                        int ii = i + filter.shift_w[x];
                        int jj = j + filter.shift_h[y];
                        if (image.valid(ii, jj)) {
                            out(i, j, b) += image(ii, jj, b) * filter.feat[b_row + f];
                        }
                        f++;
                    }
                }
            }
        }
    }

    return out;
}


Filter Filter::gaussian(int size, float sigma)
{
    Filter filter(size, size, 1);

    const float var = 2.0f * sigma * sigma;

    float sum = 0;
    for (int j = 0; j < filter.h; j++) {
        for (int i = 0; i < filter.w; i++) {
            int x = filter.shift_w[i];
            int y = filter.shift_h[j];
            float w = expf(-(x * x + y * y) / var);
            filter(i, j) = w;
            sum += w;
        }
    }

    for (int i = 0; i < filter.w * filter.h; i++) {
        filter.feat[i] /= sum;
    }

    return filter;
}


Image<float> Filter::anisotropic(const Image<float> &image, conduction_function fun_type,
                                 int iters, float lambda, float kappa)
{
    Image<float> out(image);
    Image<float> update(image);

    double (*cond_fun)(const float *, const float *, int, float) = nullptr;

    switch (fun_type) {
        case exponential:
            cond_fun = &Filter::expConduc;
            break;
        case quadratic:
            cond_fun = &Filter::quadConduct;
            break;
    }

    #ifdef _OPENMP
        #pragma omp parallel
    #endif
    for (int it = 0; it < iters; it++)
    {
        update.setFeats(out.getFeats());
        #ifdef _OPENMP
            #pragma omp for
        #endif
        for (int j = 0; j < out.getHeight(); j++) {
            for (int i = 0; i < out.getWidth(); i++)
            {
                double coeff[4] = { };

                if (out.valid(i + 1, j)) {
                    coeff[0] = cond_fun(update.getFeats(i, j), update.getFeats(i + 1, j), out.getBands(), kappa);
                }

                if (out.valid(i, j + 1)) {
                    coeff[1] = cond_fun(update.getFeats(i, j), update.getFeats(i, j + 1), out.getBands(), kappa);
                }

                if (out.valid(i - 1, j)) {
                    coeff[2] = cond_fun(update.getFeats(i, j), update.getFeats(i - 1, j), out.getBands(), kappa);
                }

                if (out.valid(i, j - 1)) {
                    coeff[3] = cond_fun(update.getFeats(i, j), update.getFeats(i, j - 1), out.getBands(), kappa);
                }

                for (int b = 0; b < out.getBands(); b++)
                {
                    if (coeff[0] > 0.0)
                        out(i, j, b) += lambda * coeff[0] * (update(i + 1, j, b) - update(i, j, b));

                    if (coeff[1] > 0.0)
                        out(i, j, b) += lambda * coeff[1] * (update(i, j + 1, b) - update(i, j, b));

                    if (coeff[2] > 0.0)
                        out(i, j, b) += lambda * coeff[2] * (update(i - 1, j, b) - update(i, j, b));

                    if (coeff[3] > 0.0)
                        out(i, j, b) += lambda * coeff[3] * (update(i, j - 1, b) - update(i, j, b));
                }
            }
        }
    }

    return out;
}