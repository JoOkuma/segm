
#ifndef SEGM_CLOSESTDYNTREE_H
#define SEGM_CLOSESTDYNTREE_H

#include "segmentation/dyntree.h"
#include <vector>

namespace segm {

    class ClosestDynTree : public DynTree
    {
    public:
        ClosestDynTree(int width, int height, int bands);
        ClosestDynTree(int width, int height, int bands, const float *feats);
        explicit ClosestDynTree(const Image<float> &image);

    protected:
        void init(const Image<int> &markers) override;
        void conquer(int x, int y, int adj_x, int adj_y) override;

    private:
        std::vector<std::vector<int>> actives; // active roots
        std::vector<std::vector<float>> dists; // active distances
    };

}


#endif //SEGM_CLOSESTDYNTREE_H
