
#include "optimumpathforest.h"

using namespace segm;

OptimumPathForest::OptimumPathForest(int n_nodes, int dim) :
        Graph(n_nodes, dim),
        heap(n_nodes),
        cost(static_cast<unsigned long>(n_nodes)),
        root(static_cast<unsigned long>(n_nodes)),
        pred(static_cast<unsigned long>(n_nodes)),
        label(static_cast<unsigned long>(n_nodes)),
        order(static_cast<unsigned long>(n_nodes))
{
    heap.setValues(cost.data());
}


void OptimumPathForest::addRoot(int node, int _label, float initial_cost)
{
    cost[node] = initial_cost;
    root[node] = node;
    label[node] = _label;
    pred[node] = nil;
    order[node] = nil;
    heap.insert(node);
}


void OptimumPathForest::run()
{
    int count = 0;
    while (!heap.isEmpty())
    {
        int p = heap.pop();
        Node &node = nodes[p];

        order[p] = count;
        count++;

        if (root[p] == p) {
            cost[p] = 0;
        }

        for (Arc &arc : node.adj)
        {
            float tmp = ((cost[p] > arc.weight) ? cost[p] : arc.weight);
            if (tmp < cost[arc.node])
            {
                cost[arc.node] = tmp;
                root[arc.node] = root[p];
                pred[arc.node] = p;
                label[arc.node] = label[p];

                if (heap.is(arc.node, heap.gray))
                    heap.goUp(arc.node);
                else
                    heap.insert(arc.node);
            }
        }
    }

    executed = true;
}


void OptimumPathForest::reset()
{
    heap.reset();

    for (int i = 0; i < size; i++)
    {
        cost[i] = 0;
        root[i] = nil;
        pred[i] = nil;
        root[i] = nil;
        label[i] = nil;
        order[i] = nil;
    }
}