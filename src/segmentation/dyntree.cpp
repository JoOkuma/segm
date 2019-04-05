
#include "dyntree.h"

using namespace segm;

DynTree::DynTree(int width, int height, int bands) :
        ForestingTransform(width, height, bands)
{
    sets = new DynSet * [width * height]();
}


DynTree::DynTree(int width, int height, int bands, const float *feats) :
        DynTree(width, height, bands)
{
    setFeats(feats);
}


DynTree::DynTree(const Image<float> &image) :
        ForestingTransform(image)
{
    sets = new DynSet * [image.getWidth() * image.getHeight()]();
}


DynTree::~DynTree()
{
    for (int i = 0; i < h * w; i++) {
        if (sets[i])
            delete sets[i];
    }
    delete[] sets;
}


DynTree::DynSet::DynSet(int bands)
{
    size = 0;
    means = new double[bands]();
}


DynTree::DynSet::~DynSet()
{
    delete[] means;
}


Image<float> DynTree::getMeans() const
{
    if (!executed)
        throw std::runtime_error("IFT must be executed before getting means");

    Image<float> output(w, h, b);
    for (int i = 0; i < w * h; i++) {
        if (sets[i]) {
            float *f = output.getFeats(i);
            for (int j = 0; j < b; j++) {
                f[j] = (float) sets[i]->means[j];
            }
        }
    }

    return output;
}


void DynTree::reset()
{
    ForestingTransform::reset();

    for (int i = 0; i < w * h; i++) {
        if (sets[i]) {
            delete sets[i];
            sets[i] = nullptr;
        }
    }
}


void DynTree::init(const Image<int> &markers)
{
    for (int p = 0; p < w * h; p++) {
        if (markers(p) >= 0) {
            cost(p) = 0;
            root(p) = p;
            label(p) = markers(p);
            heap.insert(p);
            sets[p] = new DynSet(b);
        }
    }
}


void DynTree::updatePath(int current)
{
    insert(sets[root(current)], current);
}


void DynTree::conquer(int x, int y, int adj_x, int adj_y)
{
    if (valid(adj_x, adj_y) && !heap.is(adj_x, adj_y, heap.black))
    {
        int p = index(x, y);
        int q = index(adj_x, adj_y);
        auto arc_weight = static_cast<float>( squaredDist(sets[root(p)], q) );
        float tmp = max(arc_weight, cost(p));
        if (tmp < cost(q))
        {
            cost(q)  = tmp;
            root(q)  = root(p);
            pred(q)  = p;
            label(q) = label(p);

            if (heap.is(q, heap.gray))
                heap.goUp(adj_x, adj_y);
            else
                heap.insert(q);
        }
    }
}


void DynTree::insert(segm::DynTree::DynSet *set, int p)
{
    set->size++;
    float *f = getFeats(p);
    for (int i = 0; i < b; i++)
        set->means[i] += (f[i] - set->means[i]) / set->size;
}


double DynTree::squaredDist(const segm::DynTree::DynSet *set, int p) const
{
    double dist = 0;
    float *f = getFeats(p);
    for (int i = 0; i < b; i++) {
        double diff = set->means[i] - f[i];
        dist += diff * diff;
    }

    return dist;
}


