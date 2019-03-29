
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
        label(width, height),
        order(width, height)
{
    heap.setValues(cost.getFeats());
}


ForestingTransform::ForestingTransform(int width, int height, int bands, const float *feats) :
        ForestingTransform(width, height, bands) {
    setFeats(feats);
}


void ForestingTransform::run(Image<int> &markers)
{
    if (markers.getWidth() != w || markers.getHeight() != h)
        throw std::invalid_argument("Marker image and original must have same dimensions");

    reset();

    init(markers);

    int count = 0;
    while(!heap.isEmpty())
    {
        int p = heap.pop();

        order(p) = count;
        count++;

        if (root(p) == p) {
            // TODO
            //  - verificar se isto não buga o heap (acredito que não, ele já saiu da fila)
            cost(p) = 0;
        }

        updatePath(p);

        Pixel pix = coord(p);
        // adjacency with radius 1
        conquer(pix.x, pix.y, pix.x + 1, pix.y);
        conquer(pix.x, pix.y, pix.x, pix.y + 1);
        conquer(pix.x, pix.y, pix.x - 1, pix.y);
        conquer(pix.x, pix.y, pix.x, pix.y - 1);
    }

    executed = true;
}


Image<float> ForestingTransform::getCost() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting costs");
    return cost.copy();
}


Image<int> ForestingTransform::getRoot() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting roots");
    return root.copy();
}


Image<int> ForestingTransform::getPred() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting predecessors");
    return pred.copy();
}


Image<int> ForestingTransform::getLabel() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting labels");
    return label.copy();
}


Image<int> ForestingTransform::getOrder() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting order");
    return order.copy();
}


Image<int> ForestingTransform::getBorder() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting border");

    Image<int> border(w, h);

    for (int i = 1; i < (w - 1); i++) {
        for (int j = 1; j < (h - 1); j++)
        {
            int lb = label(i, j);
            if (lb != label(i + 1, j + 1) || lb != label(i - 1, j + 1) ||
                lb != label(i - 1, j - 1) || lb != label(i + 1, j - 1) ||
                lb != label(i + 1, j) || lb != label(i, j + 1) ||
                lb != label(i - 1, j) || lb != label(i, j - 1))
            {
                border(i, j) = 1;
            }
        }
    }

    return border;
}


Image<int> ForestingTransform::getPredCount() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting predecessors count");

    // TODO
    //  - faster implementation (worst case of this is O(n^2)
    Image<int> count(w, h);
    for (int p = 0; p < w * h; p++) {
        for (int q = p; pred(q) != nil; q = pred(q)) {
            count(pred(q)) += 1;
        }
    }

    return count;
}


Image<int> ForestingTransform::getLeafPredCount() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting predecessors count");

    Image<bool> control(w, h);

    for (int p = 0; p < w * h; p++) {
        control(p) = true;
    }

    // selecting leafs
    for (int p = 0; p < w * h; p++) {
        if (pred(p) != nil) {
            control(pred(p)) = false;
        }
    }

    // TODO
    //  - faster implementation (worst case of this is O(n^2)
    Image<int> count(w, h);
    for (int p = 0; p < w * h; p++) {
        if (control(p)) {
            for (int q = p; pred(q) != nil; q = pred(q)) {
                count(pred(q)) += 1;
            }
        }
    }

    return count;
}


void ForestingTransform::trim(int index)
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before pruning trees");

    auto *float_order = new float[h * w];
    auto *color = new int[h * w]();
    Heap heap(w, h, float_order);

    const int r = root(index);
    const int o = order(index);
    color[index] = 1; // marking the correct tree

    for (int p = 0; p < h * w; p++) {
        if (root(p) == r && order(p) > o) {
            float_order[p] = order(p);
            heap.insert(p);
        }
    }

    while (!heap.isEmpty())
    {
        int p = heap.pop();

        if (pred(p) == nil)
            continue;

        if (color[pred(p)])
        {
            cost(p)  = std::numeric_limits<float>::max();
            pred(p)  = nil;
            root(p)  = nil;
            label(p) = nil;
            order(p) = nil;
            color[p] = 1;
        }
    }

    delete[] float_order;
    delete[] color;
}


void ForestingTransform::trim(int x, int y)
{
    trim(index(x, y));
}


void ForestingTransform::reset()
{
    heap.reset();

    for (int p = 0; p < w * h; p++) {
        cost(p)  = std::numeric_limits<float>::max();
        pred(p)  = nil;
        root(p)  = nil;
        label(p) = nil;
        order(p) = nil;
    }
}


void ForestingTransform::init(const Image<int> &markers)
{
    for (int p = 0; p < w * h; p++) {
        if (markers(p) >= 0) {
            cost(p) = 0;
            root(p) = p;
            label(p) = markers(p);
            heap.insert(p);
        }
    }
}
