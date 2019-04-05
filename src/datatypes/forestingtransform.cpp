
#include <stdexcept>
#include <limits>

#include "forestingtransform.h"
#include "utils/algorithm.h"

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


ForestingTransform &ForestingTransform::operator=(const segm::ForestingTransform &forest)
{
    if (w != forest.getWidth() || h != forest.getHeight() || b != forest.getBands())
    {
        if (allocated)
            delete[] feat;
        delete[] row;
        delete[] pos;
        delete[] row_index;

        w = forest.getWidth();
        h = forest.getHeight();
        b = forest.getBands();
        allocated = true;

        feat = new float[w * h * b];
        row = new int[h];
        row_index = new int[h];
        pos = new int[w * h];
        for (int i = 0, r = 0; i < h; i++, r += w) {
            row[i] = r * b;
            row_index[i] = r;
        }

        for (int i = 0, c = 0; i < w * h; i++, c += b)
            pos[i] = c;

        heap = Heap(w, h);
        cost = Image<float>(w, h);
        root = Image<int>(w, h);
        pred = Image<int>(w, h);
        label = Image<int>(w, h);
        order = Image<int>(w, h);

        executed = false;
    }

    setFeats(forest.getFeats());
    heap.setValues(cost.getFeats());

    return (*this);
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

    for (int y = 1; y < (h - 1); y++) {
        for (int x = 1; x < (w - 1); x++) {
            int lb = label(x, y);
            if (lb != label(x + 1, y + 1) || lb != label(x - 1, y + 1) ||
                lb != label(x - 1, y - 1) || lb != label(x + 1, y - 1) ||
                lb != label(x + 1, y) || lb != label(x, y + 1) ||
                lb != label(x - 1, y) || lb != label(x, y - 1)) {
                border(x, y) = 1;
            }
        }
    }

    return border;
}


Image<int> ForestingTransform::getPredCount() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting predecessors count");

    std::vector<int> order_vec(static_cast<size_t >(w * h));
    order_vec.assign(order.getFeats(), order.getFeats() + w * h);
    std::vector<int> index = indexesDecreasing(order_vec);

    Image<int> count(w, h);
    for (int i = 0; i < w * h; i++) {
        int p = index[i];
        if (pred(p) != nil)
            count(pred(p)) += count(p) + 1;
    }

    return count;
}


Image<int> ForestingTransform::getLeafPredCount() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting predecessors count");

    Image<bool> leafs(w, h);

    for (int p = 0; p < w * h; p++) {
        leafs(p) = true;
    }

    // selecting leafs
    for (int p = 0; p < w * h; p++) {
        if (pred(p) != nil) {
            leafs(pred(p)) = false;
        }
    }

    std::vector<int> order_vec(static_cast<size_t >(w * h));
    order_vec.assign(order.getFeats(), order.getFeats() + w * h);
    std::vector<int> index = indexesDecreasing(order_vec);

    Image<int> count(w, h);
    for (int i = 0; i < w * h; i++) {
        int p = index[i];
        if (pred(p) != nil)
            count(pred(p)) += count(p) + leafs(p);
    }

    return count;
}


Image<bool> ForestingTransform::getBranch(int index)
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting tree branch.");

    std::vector<int> order_vec(static_cast<size_t >(w * h));
    order_vec.assign(order.getFeats(), order.getFeats() + w * h);
    std::vector<int> index_vec = indexesIncreasing(order_vec);

    Image<bool> color(w, h);
    color(index) = true;

    for (int i = 0; i < h * w; i++) {
        int p = index_vec[i];
        if (pred(p) != nil && color(pred(p)))
            color(p) = true;
    }
    
    color(index) = false;
    return color;
}


Image<bool> ForestingTransform::getBranch(int x, int y)
{
    return getBranch(index(x, y));
}


void ForestingTransform::trim(int index)
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before pruning trees");

    std::vector<int> order_vec(static_cast<size_t >(w * h));
    order_vec.assign(order.getFeats(), order.getFeats() + w * h);
    std::vector<int> index_vec = indexesIncreasing(order_vec);

    Image<bool> color(w, h);
    color(index) = true;

    for (int i = 0; i < h * w; i++) {
        int p = index_vec[i];
        if (pred(p) != nil && color(pred(p)))
        {
            cost(p)  = 0;
            pred(p)  = nil;
            root(p)  = nil;
            label(p) = nil;
            order(p) = nil;
            color(p) = true;
        }
    }
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
