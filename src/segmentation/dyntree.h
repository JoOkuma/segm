#ifndef SEGM_DYNTREE_H
#define SEGM_DYNTREE_H

#include "datatypes/forestingtransform.h"

namespace segm
{
    class DynTree : public virtual ForestingTransform
    {

    public:
        DynTree(int width, int height, int bands);
        DynTree(int width, int height, int bands, const float *feats);
        explicit DynTree(const Image<float> &image);
        ~DynTree() override;

    protected:
        struct DynSet {
            double *means;
            int size;

            explicit DynSet(int bands);
            ~DynSet();
        };

        DynSet **sets;

    public:
        Image<float> getMeans() const;

    protected:
        void reset() override;
        void init(const Image<int> &markers) override;
        void updatePath(int current) override;

        double squaredDist(const DynSet *set, int p) const;

    private:
        void conquer(int x, int y, int adj_x, int adj_y) override;
        void insert(DynSet *set, int p);
    };

}


#endif //SEGM_DYNTREE_H
