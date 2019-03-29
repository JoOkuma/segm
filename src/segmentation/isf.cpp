
#include "datatypes/heap.h"
#include "segmentation/isf.h"

#include <eigen3/Eigen/Dense>

#include <cmath>
#include <queue>

using namespace segm;


void ISF::run(int n_superpixels, int iterations)
{
    Image<int> markers(w, h);
    for (int ite = 0; ite < iterations; ite++)
    {
        if (ite == 0)
            markers = sample(n_superpixels);
        else
            markers = computeCentroids();

        ForestingTransform::run(markers);
    }
}


void ISF::conquer(int x, int y, int adj_x, int adj_y)
{
    if (valid(adj_x, adj_y) && !heap.is(adj_x, adj_y, heap.black))
    {
        int p = index(x, y);
        int q = index(adj_x, adj_y);
        float arc_weight = cost(p) + powf(alpha * l2norm(root(p), q), beta) + 1.0f;
        if (arc_weight < cost(q))
        {
            cost(q)  = arc_weight;
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


Image<int> ISF::sample(int sample_size)
{
    auto step_x = static_cast<int>(w / round(sqrt(sample_size)));
    auto step_y = static_cast<int>(h / round(sqrt(sample_size)));

    if (step_x < 1.0f || step_y < 1.0f)
        std::runtime_error("Sample size is too big., ISF::gridSample");

    Image<int> samples(w, h);
    samples.fill(-1);

    int label = 0;
    for (int x = step_x / 2; x < w; x += step_x) {
        for (int y = step_y / 2; y < h; y += step_y) {
            samples(x, y) = label;
            label++;
        }
    }

    return samples;
}


Image<int> ISF::computeCentroids()
{
    int n_sup = label.max() + 1;

    /* 0 = x-axis, 1 = y-axis, 2 = count */
    Eigen::MatrixXi centroids = Eigen::MatrixXi::Constant(n_sup, 3, 0);

    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            int lb = label(x, y);
            centroids(lb, 0) += x;
            centroids(lb, 1) += y;
            centroids(lb, 2) += 1;
        }
    }

    Image<int> samples(w, h);
    samples.fill(-1);
    for (int lb = 0; lb < centroids.rows(); lb++)
    {
        int x = centroids(lb, 0) /= centroids(lb, 2);
        int y = centroids(lb, 1) /= centroids(lb, 2);
        if (label(x, y) == lb)
            samples(x, y) = lb;
        else {
            Pixel p = findNearest(x, y, lb);
            samples(p.x, p.y) = lb;
        }
    }

    return samples;
}


ISF::Pixel ISF::findNearest(int x, int y, int _label)
{
    Pixel p(x, y);
    std::queue<Pixel> Q;
    while (label(p.x, p.y) != _label)
    {
        if (valid(p.x + 1, p.y))
            Q.push(Pixel(p.x + 1, p.y));

        if (valid(p.x, p.y + 1))
            Q.push(Pixel(p.x, p.y + 1 ));

        if (valid(p.x - 1, p.y))
            Q.push(Pixel(p.x - 1, p.y));

        if (valid(p.x, p.y - 1))
            Q.push(Pixel(p.x, p.y - 1));

        p = Q.front();
        Q.pop();
    }

    return p;
}