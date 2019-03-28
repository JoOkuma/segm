
#ifndef SEGM_OPTIMUMPATHFOREST_H
#define SEGM_OPTIMUMPATHFOREST_H

#include <datatypes/heap.h>
#include <datatypes/graph.h>

namespace segm {

    class OptimumPathForest : public Graph<float>
    {

    public:
        OptimumPathForest(int n_nodes, int dim);

        void addRoot(int node, int label, float initial_cost = 0);

        void run();
        void reset();

    public:

        static const int nil = -1;

    protected:

        Heap heap;
        std::vector<float> cost;
        std::vector<int> root;
        std::vector<int> pred;
        std::vector<int> label;
        std::vector<int> order;

        bool executed;
    };

}


#endif //SEGM_OPTIMUMPATHFOREST_H
