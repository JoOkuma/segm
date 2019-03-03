
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

        inline float squaredl2norm(int p, int q) const
        {
            float dist = 0;
            float *v1 = getFeats(p);
            float *v2 = getFeats(q);

            for (int i = 0; i < b; i++) {
                float diff = v1[i] - v2[i];
                dist += diff * diff;
            }

            return dist;
        }
    };
}


#endif //SEGM_WATERCUT_H
