
#ifndef SEGM_SEGMIMAGE_H
#define SEGM_SEGMIMAGE_H

#include <stdexcept>
#include <limits>

#include "utils/color.h"

namespace segm {
    template<typename T>
    class Image {

    public:

        typedef struct {
            int x;
            int y;
        } Pixel;

        typedef enum : unsigned int {
            rgb = 0x00,
            xyz = 0x01,
            lab = 0x02,
            ypbpr = 0x03,
            hsv = 0x04,
            gray = 0x05,
            rgchroma = 0x10
        } ColorSpace;

        Image(int width, int height);
        Image(int width, int height, int bands);
        Image(int width, int height, int bands, const T *_feat);
        Image(int width, int height, const T *_feat);
        Image(int width, int height, int bands, T *_feat, bool alloc = true);
        Image(int width, int height, T *_feat, bool alloc = true);
        Image(const Image<T> &image);

        virtual ~Image();

        int getHeight() const { return h; }

        int getWidth() const { return w; }

        int getSize() const { return w * h; }

        int getBands() const { return b; }

        T *getFeats(int x, int y) const { return &feat[row[y] + pos[x]]; }

        T *getFeats(int p) const { return &feat[pos[p]]; }

        T *getFeats() const { return feat; }

        void setFeats(const T *array);

        Pixel coord(int p) const {
            return Pixel{.x = p % w, .y = p / w};
        }

        int index(int x, int y) const { return row_index[y] + x; }

        T max(T a, T b) { return ((a > b) ? a : b); }

        Image<T> &operator=(const Image<T> &image);

        T &operator()(int x, int y, int b) { return feat[row[y] + pos[x] + b]; }
        T operator()(int x, int y, int b) const { return feat[row[y] + pos[x] + b]; }

        /* only for 1 band images */
        T &operator()(int x, int y) { return feat[row[y] + x]; }
        T operator()(int x, int y) const { return feat[row[y] + x]; }

        T &operator()(int p) { return feat[p]; }
        T operator()(int p) const { return feat[p]; }

        inline bool valid(int x, int y) const { return ((x >= 0 && x < w) && (y >= 0 && y < h)); }

        Image<T> copy() const;
		void fill(T value);

        T max() const;
        T min() const;

        int argmax() const;
        int argmin() const;

        template<typename U>
        Image<U> convert() const;

        Image<double> convert(ColorSpace from, ColorSpace to, double normalization = 1) const;

        template<typename U>
        Image<U> rescale() const;

        bool isChromatic(ColorSpace space) { return (space & 0x10); }

    protected:
        int h = 0; /* h */
        int w = 0; /* w */
        int b = 0; /* bands */
        T *feat = nullptr;

        /* lookup table */
        int *row = nullptr; /* w plus band length padding */
        int *pos = nullptr; /* band padding only */
        int *row_index = nullptr; /* w padding only */

        bool allocated = false; /* option to not alloc features, just point to it */
    };

    template<typename T>
    Image<T>::Image(int width, int height, int bands) {
        w = width;
        h = height;
        b = bands;
        allocated = true;

        feat = new T[w * h * b]();
        row = new int[h];
        row_index = new int[h];
        pos = new int[w * h];
        for (int i = 0, r = 0; i < h; i++, r += w) {
            row[i] = r * b;
            row_index[i] = r;
        }
        for (int i = 0, c = 0; i < w * h; i++, c += b)
            pos[i] = c;
    }

    template<typename T>
    Image<T>::Image(int width, int height, int bands, const T *_feat) {
        w = width;
        h = height;
        b = bands;
        allocated = true;

        feat = new T[w * h * b];
        row = new int[h];
        row_index = new int[h];
        pos = new int[w * h];
        for (int i = 0, r = 0; i < h; i++, r += w) {
            row[i] = r * b;
            row_index[i] = r;
        }

        for (int i = 0, c = 0; i < w * h; i++, c += b)
            pos[i] = c;

        setFeats(_feat);
    }

    template<typename T>
    Image<T>::Image(int width, int height, int bands, T *_feat, bool alloc) {
        w = width;
        h = height;
        b = bands;

        if (alloc) {
            allocated = true;
            feat = new T[w * h * b];
            setFeats(_feat);
        } else {
            allocated = false;
            feat = _feat;
        }

        row = new int[h];
        row_index = new int[h];
        pos = new int[w * h];
        for (int i = 0, r = 0; i < h; i++, r += w) {
            row[i] = r * b;
            row_index[i] = r;
        }

        for (int i = 0, c = 0; i < w * h; i++, c += b)
            pos[i] = c;
    }

    template<typename T>
    Image<T>::Image(int width, int height) :
            Image<T>(width, height, 1) { }

    template<typename T>
    Image<T>::Image(int width, int height, const T *_feat) :
            Image<T>(width, height, 1, _feat) { }

    template<typename T>
    Image<T>::Image(int width, int height, T *_feat, bool alloc) :
            Image<T>(width, height, 1, _feat, alloc) { }

    template<typename T>
    Image<T>::Image(const Image<T> &image) :
            Image<T>(image.getWidth(), image.getHeight(), image.getBands(), image.getFeats()) { }

    template<typename T>
    Image<T>::~Image() {
        if (allocated)
            delete[] feat;
        delete[] row;
        delete[] pos;
        delete[] row_index;
    };

    template<typename T>
    void Image<T>::setFeats(const T *array)
    {
        for (int i = 0; i < w * h * b; i++)
            feat[i] = array[i];
    }

    template<typename T>
    Image<T> &Image<T>::operator=(const Image<T> &image)
    {
        if (w != image.getWidth() || h != image.getHeight() ||
            b != image.getBands() || !allocated)
        {
            if (allocated)
                delete[] feat;
            delete[] row;
            delete[] pos;
            delete[] row_index;

            w = image.getWidth();
            h = image.getHeight();
            b = image.getBands();
            allocated = true;

            feat = new T[w * h * b];
            row = new int[h];
            row_index = new int[h];
            pos = new int[w * h];
            for (int i = 0, r = 0; i < h; i++, r += w) {
                row[i] = r * b;
                row_index[i] = r;
            }

            for (int i = 0, c = 0; i < w * h; i++, c += b)
                pos[i] = c;
        }

        setFeats(image.getFeats());

        return (*this);
    }

    template<typename T>
    Image<T> Image<T>::copy() const {
        Image<T> out(w, h, b);
        out.setFeats(feat);
        return out;
    }

	template<typename T>
	void Image<T>::fill(T value) {
		for (int i = 0; i < w * h * b; i++)
			feat[i] = value;
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


    template<typename T>
    template<typename U>
    Image<U> Image<T>::convert() const
    {
        Image<U> out(w, h, b);
        for (int p = 0; p < w * h * b; p++) {
            out(p) = static_cast<U>(feat[p]);
        }

        return out;
    }


    template<typename T>
    Image<double> Image<T>::convert(Image<T>::ColorSpace from, Image<T>::ColorSpace to, double normalization) const
    {
        if (from == to)
            return convert<double>();

        void (*convFun)(const double *, double*) = nullptr;

        unsigned int conversion = (from << 16) | to;
        switch (conversion)
        {
            case (rgb << 16) | xyz:
                convFun = rgb2xyz;
                break;
            case (xyz << 16) | rgb:
                convFun = xyz2rgb;
                break;
            case (xyz << 16) | lab:
                convFun = xyz2lab;
                break;
            case (lab << 16) | xyz:
                convFun = lab2xyz;
                break;
            case (rgb << 16) | lab:
                convFun = rgb2lab;
                break;
            case (lab << 16) | rgb:
                convFun = lab2rgb;
                break;
            case (rgb << 16) | ypbpr:
                convFun = rgb2ypbpr;
                break;
            case (ypbpr << 16) | rgb:
                convFun = ypbpr2rgb;
                break;
            case (rgb << 16) | hsv:
                convFun = rgb2hsv;
                break;
            case (hsv << 16) | rgb:
                convFun = hsv2rgb;
                break;
            case (rgb << 16) | gray:
                convFun = rgb2gray;
                break;
            case (gray << 16) | rgb:
                convFun = gray2rgb;
                break;
            case (rgb << 16) | rgchroma:
                convFun = rgb2rgchroma;
                break;
            default:
                throw std::invalid_argument("Color conversion requested not found");
        }

        double dbl_feat[3];
        Image<double> out(w, h, ((to != gray) ? 3 : 1));

        for (int i = 0, p = 0; i < w * h; i++) {
            switch (from) {
                case gray:
                    dbl_feat[0] = feat[p] / normalization;
                    convFun(dbl_feat, out.getFeats(p));
                    p++;
                    break;
                default:
                    dbl_feat[0] = feat[p] / normalization;
                    dbl_feat[1] = feat[p + 1] / normalization;
                    dbl_feat[2] = feat[p + 2] / normalization;
                    convFun(dbl_feat, out.getFeats(i));
                    p += 3;
                    break;
            }
        }

        return out;
    }


    template<typename T>
    template<typename U>
    Image<U> Image<T>::rescale() const
    {
        Image<U> out(w, h, b);
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for (int bb = 0; bb < b; bb++) {
            T max = std::numeric_limits<T>::min();
            T min = std::numeric_limits<T>::max();
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    T val = (*this)(i, j, bb);
                    if (val > max) max = val;
                    if (val < min) min = val;
                }
            }

            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    out(i, j, bb) = (U) (((*this)(i, j, bb) - min) / max);
                }
            }
        }
        return out;
    }

}
#endif //SEGM_SEGMIMAGE_H

