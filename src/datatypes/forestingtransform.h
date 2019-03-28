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
        explicit ForestingTransform(const Image<float> &image) :
                ForestingTransform(image.getWidth(), image.getHeight(), image.getBands(), image.getFeats()) { }
        ~ForestingTransform() override = default;

        void run(Image<int> &markers);

        Image<float> getCost() const;
        Image<int> getRoot() const;
        Image<int> getPred() const;
        Image<int> getLabel() const;
        Image<int> getOrder() const;

        Image<int> getPredCount() const;
        Image<int> getLeafPredCount() const;

        void trim(int index);
        void trim(int x, int y);

    public:

        static const int nil = -1;

    protected:

        virtual void reset();
        virtual void init(const Image<int> &markers);
        virtual void updatePath(int) { };
        virtual void conquer(int x, int y, int adj_x, int adj_y) = 0;

    protected:

        Heap heap;
        Image<float> cost;
        Image<int> root;
        Image<int> pred;
        Image<int> label;
        Image<int> order;

        bool executed = false;

    };
}


#endif //SEGM_FORESTINGTRANSFORM_H
