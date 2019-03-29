
#ifndef SEGM_ISF_H
#define SEGM_ISF_H

#include "datatypes/forestingtransform.h"

#include <cmath>

namespace segm {

    class ISF : public ForestingTransform 
    {

    public:
        ISF(int width, int height, int bands, float _alpha = 5.0f, float _beta = 12.0f) :
                ForestingTransform(width, height, bands), alpha(_alpha), beta(_beta) { }
        ISF(int width, int height, int bands, const float *feats, float _alpha = 5.0f, float _beta = 12.0f) :
                ForestingTransform(width, height, bands, feats), alpha(_alpha), beta(_beta) { }
        explicit ISF(const Image<float> &image, float _alpha = 5.0f, float _beta = 12.0f) :
                ForestingTransform(image), alpha(_alpha), beta(_beta) { }

        void run(int n_superpixels, int iterations = 5);

    private:

        void conquer(int x, int y, int adj_x, int adj_y) override;

        inline float l2norm(int p, int q) const
        {
            float dist = 0;
            float *v1 = getFeats(p);
            float *v2 = getFeats(q);

            for (int i = 0; i < b; i++) {
                float diff = v1[i] - v2[i];
                dist += diff * diff;
            }

            return sqrtf(dist);
        }

        virtual Image<int> sample(int sample_size); /* default is grid sample */
        virtual Image<int> computeCentroids(); /* default is geodesic center */
        Pixel findNearest(int x, int y, int _label);

    private:

        float alpha;
        float beta;

    };

}


#endif //SEGM_ISF_H
