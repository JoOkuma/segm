
#include "closestdyntree.h"

using namespace segm;


ClosestDynTree::ClosestDynTree(int width, int height, int bands) :
        DynTree(width, height, bands) { }


ClosestDynTree::ClosestDynTree(int width, int height, int bands, const float *feats) :
        DynTree(width, height, bands, feats) { }


ClosestDynTree::ClosestDynTree(const Image<float> &image) :
        DynTree(image) { }


void ClosestDynTree::init(const Image<int> &markers)
{
    size_t n_labels = (size_t) markers.max() + 1;
    actives.resize(n_labels);
    dists.resize(n_labels);
    for (int p = 0; p < w * h; p++) {
        if (markers(p) >= 0) {
            cost(p) = 0;
            root(p) = p;
            label(p) = markers(p);
            heap.insert(p);
            sets[p] = new DynSet(b);
            actives.at((size_t) markers(p)).push_back(p);
            dists.at((size_t) markers(p)).push_back(0);
        }
    }
}


void ClosestDynTree::conquer(int x, int y, int adj_x, int adj_y)
{
    if (valid(adj_x, adj_y) && !heap.is(adj_x, adj_y, heap.black))
    {
        int p = index(x, y);
        int q = index(adj_x, adj_y);
        float arc_weight = std::numeric_limits<float>::max();
        int min_r = -1;

        std::vector<int> &active = actives[label(p)];
        std::vector<float> &dist = dists[label(p)];
//        #ifdef _OPENMP
//            #pragma omp parallel for
//        #endif
        for (size_t i = 0; i < active.size(); i++) {
            dist[i] = (float) squaredDist(sets[active[i]], q);
        }

        for (size_t i = 0; i < active.size(); i++) {
            if (dist[i] < arc_weight) {
                min_r = active[i];
                arc_weight = dist[i];
            }
        }

        float tmp = max(arc_weight, cost(p));
        if (tmp < cost(q))
        {
            cost(q)  = tmp;
            root(q)  = min_r;
            pred(q)  = p;
            label(q) = label(p);

            if (heap.is(q, heap.gray))
                heap.goUp(adj_x, adj_y);
            else
                heap.insert(q);
        }
    }
}