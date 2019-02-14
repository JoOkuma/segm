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
        ~ForestingTransform() override;

        void run(Image<int> &markers, float plato = 0.01f);

        Image<float> getCost() const;
        Image<int> getRoot() const;
        Image<int> getPred() const;
        Image<int> getLabel() const;
        Image<int> getOrder() const;

        Image<int> getPredCount() const;
        Image<int> getLeafPredCount() const;

    protected:
        const int nil = -1;

        Heap heap;
        Image<float> cost;
        Image<int> root;
        Image<int> pred;
        Image<int> label;
        Image<int> order;

        bool executed = false;

    private:
        virtual void reset();
        virtual void updatePath() { };
        virtual void conquer(int x, int y, int adj_x, int adj_y) = 0;
    };
}


#endif //SEGM_FORESTINGTRANSFORM_H
