
#include "dtisf.h"

using namespace segm;


DTISF &DTISF::operator=(const DTISF &dtisf)
{
    if (w != dtisf.getWidth() || h != dtisf.getHeight()) {
        for (int i = 0; i < h * w; i++) {
            if (sets[i])
                delete sets[i];
        }
        delete[] sets;
        sets = new DynSet *[dtisf.getWidth() * dtisf.getWidth()]();
    }

    ISF::operator=(dtisf);

    return (*this);
}


void DTISF::run(int n_superpixels, int iterations)
{
    Image<int> markers(w, h);
    for (int ite = 0; ite < iterations; ite++)
    {
        if (ite == 0)
            markers = ISF::sample(n_superpixels);
        else
            markers = ISF::computeCentroids();

        DynTree::run(markers);
    }
}


void DTISF::conquer(int x, int y, int adj_x, int adj_y)
{
    if (valid(adj_x, adj_y) && !heap.is(adj_x, adj_y, heap.black))
    {
        int p = index(x, y);
        int q = index(adj_x, adj_y);
        auto arc_weight = cost(p) + pow(alpha * sqrt(squaredDist(sets[root(p)], q)), beta) + 1.0;
        float tmp = max(static_cast<float>(arc_weight), cost(p));
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