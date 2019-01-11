#ifndef SEGM_FORESTINGTRANSFORM_H
#define SEGM_FORESTINGTRANSFORM_H

#include "datatypes/heap.h"
#include "datatypes/image.h"

namespace segm
{
    class ForestingTransform : public Image<float>
    {

    public:
        ForestingTransform(int width, int height, int bands);
        ForestingTransform(int width, int height, int bands, const float *feats);
        ~ForestingTransform() override;

        void run(Image<int> &markers, int height = 1);

        Image<float> getCost();
        Image<int> getRoot();
        Image<int> getPred();
        Image<int> getLabel();

    protected:
        const int nil = -1;

        Heap heap;
        Image<float> cost;
        Image<int> root;
        Image<int> pred;
        Image<int> label;

        bool executed = false;

    private:
        virtual void reset();
        virtual void updatePath() { };
        virtual void conquer(int x, int y, int adj_x, int adj_y) { };
    };
}


#endif //SEGM_FORESTINGTRANSFORM_H
