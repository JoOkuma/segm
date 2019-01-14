
#include <limits>
#include <stdexcept>
#include "heap.h"

using namespace segm;

Heap::Heap(int width, int height)
{
    size = width * height;
    w = width;
    h = height;
    color = new heap_color[size];
    node = new int[size];
    pos = new int[size];
    last = -1;
    policy = min_value;
    row = new int[h];

    for (int i = 0; i < size; i++) {
        color[i] = white;
        pos[i] = -1;
        node[i] = -1;
    }

    for (int i = 0, r = 0; i < h; i++, r += w)
        row[i] = r;
}

Heap::Heap(int width, int height, float *pathval)
{
    size = width * height;
    w = width;
    h = height;
    value = pathval;
    color = new heap_color[size];
    node = new int[size];
    pos = new int[size];
    last = -1;
    policy = min_value;
    row = new int[h];

    for (int i = 0; i < size; i++) {
        color[i] = white;
        pos[i] = -1;
        node[i] = -1;
    }

    for (int i = 0, r = 0; i < h; i++, r += width)
        row[i] = r;

}

Heap::~Heap()
{
    delete[] color;
    delete[] node;
    delete[] pos;
    delete[] row;
}

void Heap::setValues(float *pathval)
{
    reset();
    value = pathval;
}

void Heap::insert(int p)
{
    if (isFull())
        throw std::runtime_error("Heap is full");

    last++;
    node[last] = p;
    color[p] = gray;
    pos[p] = last;
    goUp(last);
}

int Heap::pop()
{
    if (isEmpty())
        throw std::runtime_error("Heap is empty");

    int p = node[0];
    pos[p] = -1;
    color[p] = black;
    node[0] = node[last];
    pos[node[0]] = 0;
    node[last] = -1;
    last--;
    goDown(0);

    return p;
}

std::pair<int, int> Heap::popPair()
{
    int p = pop();
    std::pair<int, int> coord(p % w, p / h);
    return coord;
}


void Heap::remove(int x, int y)
{
    int p = row[y] + x;
    if (pos[p] == -1)
        throw std::invalid_argument("Element is not in the Heap");

    float val = value[p];

    if (policy == min_value)
        value[p] = std::numeric_limits<float>::min();
    else
        value[p] = std::numeric_limits<float>::max();

    goUp(pos[p]);
    pop();
    value[p] = val;
    color[p] = white;
}

void Heap::reset()
{
    for (int i = 0; i < size; i++) {
        color[i] = white;
        pos[i] = -1;
        node[i] = -1;
    }
    last = -1;
}

void Heap::swap(int &i, int &j)
{
    int tmp = i;
    i = j;
    j = tmp;
}

void Heap::swapUp(int &dad, int &son)
{
    swap(node[dad], node[son]);
    pos[node[son]] = son;
    pos[node[dad]] = dad;
    son = dad;
    dad = getDad(son);
}

void Heap::goUp(int position)
{
    int dad = getDad(position);
    int son = position;
    if (policy == min_value) {
        while ((dad >= 0) && (value[node[dad]] > value[node[son]])) {
            swapUp(dad, son);
        }
    } else {
        while ((dad >= 0) && (value[node[dad]] < value[node[son]])) {
            swapUp(dad, son);
        }
    }
}

void Heap::goDown(int position)
{
    int dad = position;
    int left = leftSon(dad);
    int right = rightSon(dad);
    if (policy == min_value) {
        if (left <= last && value[node[left]] < value[node[dad]])
            dad = left;
        if (right <= last && value[node[right]] < value[node[dad]])
            dad = right;
    } else {
        if (left <= last && value[node[left]] > value[node[dad]])
            dad = left;
        if (right <= last && value[node[right]] > value[node[dad]])
            dad = right;
    }

    if (position != dad) {
        swap(node[dad], node[position]);
        pos[node[dad]] = dad;
        pos[node[position]] = position;
        goDown(dad);
    }
}

