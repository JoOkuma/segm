
#ifndef SEGM_DTISF_H
#define SEGM_DTISF_H

#include "segmentation/isf.h"
#include "segmentation/dyntree.h"

namespace segm
{
    class DTISF : public DynTree, public ISF
    {

    public:
        DTISF(int width, int height, int bands, float _alpha = 5.0f, float _beta = 12.0f) :
                ForestingTransform(width, height, bands),
                DynTree(width, height, bands),
                ISF(width, height, bands, _alpha, _beta) { }
        DTISF(int width, int height, int bands, const float *feats, float _alpha = 5.0f, float _beta = 12.0f) :
                ForestingTransform(width, height, bands, feats),
                DynTree(width, height, bands, feats),
                ISF(width, height, bands, feats, _alpha, _beta) { }
        explicit DTISF(const Image<float> &image, float _alpha = 5.0f, float _beta = 12.0f) :
                ForestingTransform(image),
                DynTree(image),
                ISF(image, _alpha, _beta) { }

        DTISF &operator=(const DTISF &dtisf);

        void run(int n_superpixels, int iterations = 5);

    private:
        void conquer(int x, int y, int adj_x, int adj_y) override;

    };

}


#endif //SEGM_DTISF_H
