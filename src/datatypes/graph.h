
#ifndef SEGM_GRAPH_H
#define SEGM_GRAPH_H

#include "datatypes/image.h"

#include <vector>
#include <list>

namespace segm {

    template<typename T>
    class Graph {

    protected:

        struct Arc {
            int node;
            T weight;
            Arc() : node(-1), weight(0) { };
            Arc(int _node, T _weight) :
                    node(_node) , weight(_weight) { };
        };

        struct Node {
            std::list<Arc> adj;
        };

    public:
        explicit Graph<T>(int n_nodes, int dim);

        int getSize() const { return size; }
        int getDims() const { return dimensions; }

        void addEdge(int current, int neighbor, T weight = 0);

    protected:
        int size;
        int dimensions;

        std::vector<Node> nodes;
    };


    template<typename T>
    Graph<T>::Graph(int n_nodes, int dim) :
            nodes(static_cast<unsigned long>(n_nodes))
    {
        size = n_nodes;
        dimensions = dim;
    }

    template<typename T>
    void Graph<T>::addEdge(int current, int neighbor, T weight)
    {
        Node &node = nodes[current];
        node.adj.push_back(Arc(neighbor, weight));
    }
}


#endif //SEGM_GRAPH_H
