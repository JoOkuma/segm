//
// Created by jookuma on 10/01/19.
//

#include <stdexcept>
#include <limits>
#include "forestingtransform.h"

using namespace segm;

ForestingTransform::ForestingTransform(int width, int height, int bands) :
        Image<float>(width, height, bands),
        heap(width, height),
        cost(width, height),
        root(width, height),
        pred(width, height),
        label(width, height)
{

}


ForestingTransform::ForestingTransform(int width, int height, int bands, const float *feats) :
        Image<float>(width, height, bands, feats),
        heap(w, h),
        cost(w, h),
        root(w, h),
        pred(w, h),
        label(w, h)
{

}

ForestingTransform::~ForestingTransform() {

}

void ForestingTransform::run(Image<int> &markers, int plato)
{
    if (markers.getWidth() != w || markers.getHeight() != h)
        throw std::invalid_argument("Marker image and original must have same dimensions");

    reset();

    for (int p = 0; p < w * h; p++) {
        if (markers(p) >= 0) {
            cost(p) = plato;
            root(p) = p;
            label(p) = markers(p);
            heap.insert(p);
        }
    }

    while(!heap.isEmpty())
    {
        int p = heap.pop();

        if (root(p) == p) {
            cost(p) = 0; // TODO VERIFICAR SE ISSO NAO BUGA HEAP
        }

        updatePath();

        Pixel pix = coord(p);
        // adjacency with radius 1
        conquer(pix.x, pix.y, pix.x + 1, pix.y);
        conquer(pix.x, pix.y, pix.x, pix.y + 1);
        conquer(pix.x, pix.y, pix.x - 1, pix.y);
        conquer(pix.x, pix.y, pix.x, pix.y - 1);
    }

    executed = true;
}


Image<float> ForestingTransform::getCost()
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting costs");
    return cost;
}

Image<int> ForestingTransform::getRoot()
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting roots");
    return root;
}

Image<int> ForestingTransform::getPred()
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting predecessors");
    return pred;
}

Image<int> ForestingTransform::getLabel()
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting labels");
    return label;
}

void ForestingTransform::reset()
{
    heap.reset();

    for (int p = 0; p < w * h; p++) {
        cost(p)  = std::numeric_limits<float>::max();
        pred(p)  = nil;
        root(p)  = nil;
        label(p) = nil;
    }

    heap.setValues(cost.getFeat());
}
