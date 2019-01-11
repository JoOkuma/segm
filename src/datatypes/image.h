//
// Created by jookuma on 10/01/19.
//

#ifndef SEGM_SEGMIMAGE_H
#define SEGM_SEGMIMAGE_H

#include <stdexcept>

namespace segm {
    template<typename T>
    class Image {

    public:
        typedef struct {
            int x;
            int y;
        } Pixel;

        Image(int width, int height);
        Image(int width, int height, int bands);
        Image(int width, int height, int bands, const T *_feat);
        Image(int width, int height, const T *_feat);

        virtual ~Image();

        int getHeight() const { return h; }

        int getWidth() const { return w; }

        T *getFeat(int x, int y) const { return &feat[row[y] + col[x]]; };

        T *getFeat() const { return feat; }

        T squaredl2norm(int x1, int y1, int x2, int y2);

        Pixel coord(int p) {
            return Pixel{.x = p % w, .y = p / h};
        };

        int index(int x, int y) const { return row_index[y] + x; }

        T max(T a, T b) { return (a > b) ? a : b; }

        /* only for 1 band images */
        T &operator()(int x, int y) { return feat[row_index[y] + x]; }
        T &operator()(int p) { return feat[p]; }

        bool valid(int x, int y) const { return ((x >= 0 && x < w) && (y >= 0 && y < h)); }

    protected:
        int h; /* h */
        int w; /* w */
        int b; /* bands */
        T *feat;

        /* lookup table */
        int *row; /* w plus band length padding */
        int *col; /* band padding only */
        int *row_index; /* w padding only */
    };

    template<typename T>
    Image<T>::Image(int width, int height) {
        w = width;
        h = height;
        b = 1;

        feat = new T[w * h];
        row = new int[h];
        row_index = new int[h];
        col = new int[w];
        for (int i = 0, r = 0; i < h; i++, r += w) {
            row[i] = r;
            row_index[i] = r;
        }
        for (int i = 0; i < w; i++)
            col[i] = i;
    }

    template<typename T>
    Image<T>::Image(int width, int height, int bands) {
        w = width;
        h = height;
        b = bands;

        feat = new T[w * h * b];
        row = new int[h];
        row_index = new int[h];
        col = new int[w];
        for (int i = 0, r = 0; i < h; i++, r += w) {
            row[i] = r * b;
            row_index[i] = r;
        }
        for (int i = 0, c = 0; i < w; i++, c += b)
            col[i] = c;
    }

    template<typename T>
    Image<T>::Image(int width, int height, int bands, const T *_feat) {
        w = width;
        h = height;
        b = bands;

        feat = new T[w * h * b];
        row = new int[h];
        row_index = new int[h];
        col = new int[w];
        for (int i = 0, r = 0; i < h; i++, r += w) {
            row[i] = r * b;
            row_index[i] = r;
        }

        for (int i = 0, c = 0; i < w; i++, c += b)
            col[i] = c;

        for (int i = 0; i < w * h * b; i++)
            feat[i] = _feat[i];
    }

    template<typename T>
    Image<T>::Image(int width, int height, const T *_feat) {
        w = width;
        h = height;
        b = 1;

        feat = new T[w * h];
        row = new int[h];
        row_index = new int[h];
        col = new int[w];
        for (int i = 0, r = 0; i < h; i++, r += w) {
            row[i] = r;
            row_index[i] = r;
        }

        for (int i = 0; i < w; i++)
            col[i] = i;

        for (int i = 0; i < w * h; i++)
            feat[i] = _feat[i];
    }

    template<typename T>
    Image<T>::~Image() {
        delete feat;
        delete row;
        delete col;
        delete row_index;
    };

    template<typename T>
    T Image<T>::squaredl2norm(int x1, int y1, int x2, int y2) {
        T dist = 0;
        T *v1 = getFeat(x1, y1);
        T *v2 = getFeat(x2, y2);

        for (int i = 0; i < b; i++) {
            T diff = v1[i] - v2[i];
            dist += diff * diff;
        }

        return dist;
    }

}
#endif //SEGM_SEGMIMAGE_H

