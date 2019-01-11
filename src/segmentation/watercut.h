//
// Created by jookuma on 10/01/19.
//

#ifndef SEGM_WATERCUT_H
#define SEGM_WATERCUT_H

#include "datatypes/image.h"
#include "datatypes/heap.h"

namespace segm
{
    class WaterCut : public Image<float>
    {

    public:
        WaterCut(int width, int height, int bands);
        WaterCut(int width, int height, int bands, const float *feats);
        ~WaterCut() override;

        virtual void run(Image<int> &markers, int height = 1);

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
        void reset();
        void conquer(int x, int y, int adj_x, int adj_y);
    };
}


#endif //SEGM_WATERCUT_H
