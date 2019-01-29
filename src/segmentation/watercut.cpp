
#include "watercut.h"

using namespace segm;

void WaterCut::conquer(int x, int y, int adj_x, int adj_y)
{
    if (valid(adj_x, adj_y) && !heap.is(adj_x, adj_y, heap.black))
    {
        int p = index(x, y);
        int q = index(adj_x, adj_y);
        float arc_weight = squaredl2norm(x, y, adj_x, adj_y);
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