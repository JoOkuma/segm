
#ifndef SEGM_HEAP_H
#define SEGM_HEAP_H

#include <utility>

namespace segm
{
    class Heap
    {

    public:
        enum heap_color {
            white,
            gray,
            black
        };

        enum removal_policy {
            min_value,
            max_value
        };


        Heap(int width, int height);
        Heap(int width, int height, float *cost);
        Heap(int size);
        Heap(int size, float *cost);
        ~Heap();

        Heap &operator=(const Heap &heap);

        int getSize() const { return size; }
        int getWidth() const { return w; }
        int getHeight() const { return h; }
        float *getValues() const { return value; }

        void setValues(float *cost);

        bool isFull() const { return (last == (size - 1)); }
        bool isEmpty() const { return (last == -1); }

        void insert(int p);
        void insert(int x, int y) { insert(row[y] + x); }

        int pop();
        std::pair<int, int> popPair();

        void remove(int x, int y);

        void reset();

        bool is(int p, heap_color col) const { return (color[p] == col); }
        bool is(int x, int y, heap_color col) const {
            return (is(row[y] + x, col));
        }

        void goUp(int x, int y) { goUpPos(pos[row[y] + x]); }
        void goDown(int x, int y) { goDownPos(pos[row[y] + x]); }
        void goUp(int idx) { goUpPos(pos[idx]); }
        void goDown(int idx) { goDownPos(pos[idx]); }

        void setPolicy(removal_policy _policy) { policy = _policy; }

    private:
        inline int getDad(int pos) { return (pos - 1) / 2; }
        inline int leftSon(int pos) { return 2 * pos + 1; }
        inline int rightSon(int pos) { return 2 * pos + 2; }

        void swap(int &i, int &j);
        void swapUp(int &dad, int &son);

        void goUpPos(int position);
        void goDownPos(int position);

    private:
        int size;
        int w;
        int h;

        float *value;
        heap_color *color;
        int *node;
        int *pos;
        int last;
        removal_policy policy;

        int *row;
    };
}

#endif //SEGM_HEAP_H
