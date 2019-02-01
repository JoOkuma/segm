
#ifndef SEGM_WATERCUT_H
#define SEGM_WATERCUT_H

#include "datatypes/forestingtransform.h"

namespace segm
{
    class WaterCut : public ForestingTransform
    {

    public:
        WaterCut(int width, int height, int bands) :
                ForestingTransform(width, height, bands) { }
        WaterCut(int width, int height, int bands, const float *feats) :
                ForestingTransform(width, height, bands, feats) { }
        explicit WaterCut(const Image<float> &image) : ForestingTransform(image) { }

    private:
        void conquer(int x, int y, int adj_x, int adj_y) override;
    };
}


#endif //SEGM_WATERCUT_H
