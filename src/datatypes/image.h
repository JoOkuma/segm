
#ifndef SEGM_SEGMIMAGE_H
#define SEGM_SEGMIMAGE_H

#include <stdexcept>
#include <limits>

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

        int getBands() const { return b; }

        T *getFeats(int x, int y) const { return &feat[row[y] + col[x]]; };

        T *getFeats() const { return feat; }

        T squaredl2norm(int x1, int y1, int x2, int y2);

        Pixel coord(int p) {
            return Pixel{.x = p % w, .y = p / h};
        };

        int index(int x, int y) const { return row_index[y] + x; }

        T max(T a, T b) { return ((a > b) ? a : b); }

        T &operator()(int x, int y, int b) { return feat[row[y] + col[x] + b]; }
        T operator()(int x, int y, int b) const { return feat[row[y] + col[x] + b]; }

        /* only for 1 band images */
        T &operator()(int x, int y) { return feat[row_index[y] + x]; }
        T operator()(int x, int y) const { return feat[row_index[y] + x]; }

        T &operator()(int p) { return feat[p]; }
        T operator()(int p) const { return feat[p]; }

        bool valid(int x, int y) const { return ((x >= 0 && x < w) && (y >= 0 && y < h)); }

        Image<T> copy() const;

        T max() const;
        T min() const;

        int argmax() const;
        int argmin() const;

    protected:
        int h; /* h */
        int w; /* w */
        int b; /* bands */
        T *feat = nullptr;

        /* lookup table */
        int *row = nullptr; /* w plus band length padding */
        int *col = nullptr; /* band padding only */
        int *row_index = nullptr; /* w padding only */
    };

    template<typename T>
    Image<T>::Image(int width, int height) {
        w = width;
        h = height;
        b = 1;

        feat = new T[w * h]();
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

        feat = new T[w * h * b]();
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
        delete[] feat;
        delete[] row;
        delete[] col;
        delete[] row_index;
    };

    template<typename T>
    inline T Image<T>::squaredl2norm(int x1, int y1, int x2, int y2) {
        T dist = 0;
        T *v1 = getFeats(x1, y1);
        T *v2 = getFeats(x2, y2);

        for (int i = 0; i < b; i++) {
            T diff = v1[i] - v2[i];
            dist += diff * diff;
        }

        return dist;
    }

    template <typename T>
    Image<T> Image<T>::copy() const {
        Image<T> out(w, h);
        for (int i = 0; i < w * h * b; i++) {
            out(i) = feat[i];
        }
        return out;
    }

    template<typename T>
    T Image<T>::max() const {
        T maximum = std::numeric_limits<T>::min();
        for (int i = 0; i < w * h * b; i++) {
            if (feat[i] > maximum) maximum = feat[i];
        }
        return maximum;
    }

    template<typename T>
    T Image<T>::min() const {
        T minimum = std::numeric_limits<T>::max();
        for (int i = 0; i < w * h * b; i++) {
            if (feat[i] < minimum) minimum = feat[i];
        }
        return minimum;
    }

    template<typename T>
    int Image<T>::argmax() const {
        T maximum = std::numeric_limits<T>::min();
        int p = -1;
        for (int i = 0; i < w * h * b; i++) {
            if (feat[i] > maximum) {
                maximum = feat[i];
                p = i;
            }
        }
        return p;
    }

    template<typename T>
    int Image<T>::argmin() const {
        T minimum = std::numeric_limits<T>::max();
        int p = -1;
        for (int i = 0; i < w * h * b; i++) {
            if (feat[i] < minimum) {
                minimum = feat[i];
                p = i;
            }
        }
        return p;
    }
}
#endif //SEGM_SEGMIMAGE_H

