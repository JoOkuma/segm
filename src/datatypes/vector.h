
#ifndef SEGM_VECTOR_H
#define SEGM_VECTOR_H

#include <iostream>

namespace segm
{
    template<typename T>
    class Vector
    {
    public:
        explicit Vector(int _size);
        Vector(int _size, const T *array);
        Vector(int _size, T *array, bool alloc = true);
        ~Vector();

        T &operator[](int index) const { return feat[index]; }
        T &operator[](int index) { return feat[index]; }

        Vector<T> &operator=(const Vector<T> &vector);
        Vector<T> &operator=(T value);

        Vector<T> &operator+=(T value);
        Vector<T> &operator-=(T value);
        Vector<T> &operator*=(T value);
        Vector<T> &operator/=(T value);

        int getSize() const { return size; }

        void setFeats(const T *array);

        T * getFeats(int i = 0) const { return &feat[i]; }

        void print();

    protected:
        int size;
        bool allocated = true;

        T* feat = nullptr;

    };

    template<typename T>
    Vector<T>::Vector(int _size)
    {
        size = _size;
        feat = new T[size]();
        allocated = true;
    }

    template<typename T>
    Vector<T>::Vector(int _size, const T *array) :
        Vector(_size) {
        setFeats(array);
    }

    template<typename T>
    Vector<T>::Vector(int _size, T *array, bool alloc)
    {
        size = _size;
        if (alloc) {
            feat = new T[size]();
            allocated = true;
            setFeats(array);
        } else {
            feat = array;
        }
    }

    template<typename T>
    Vector<T>::~Vector()
    {
        if (allocated)
            delete[] feat;
    }

    template<typename T>
    Vector<T> &Vector<T>::operator=(const Vector<T> &vector)
    {
        if (size != vector.size || !allocated) {
            if (allocated)
                delete[] feat;
            size = vector.size;
            feat = new T[size]();
            allocated = true;
        }

        setFeats(vector.getFeats());

        return (*this);
    }

    template<typename T>
    Vector<T> &Vector<T>::operator=(T value) {
        for (int i = 0; i < size; i++)
            feat[i] = value;

        return (*this);
    }

    template<typename T>
    Vector<T> &Vector<T>::operator+=(T value) {
        for (int i = 0; i < size; i++)
            feat[i] += value;
        return (*this);
    }

    template<typename T>
    Vector<T> &Vector<T>::operator-=(T value) {
        for (int i = 0; i < size; i++)
            feat[i] -= value;
        return (*this);
    }

    template<typename T>
    Vector<T> &Vector<T>::operator*=(T value) {
        for (int i = 0; i < size; i++)
            feat[i] *= value;
        return (*this);
    }

    template<typename T>
    Vector<T> &Vector<T>::operator/=(T value) {
        for (int i = 0; i < size; i++)
            feat[i] /= value;
        return (*this);
    }

    template<typename T>
    void Vector<T>::setFeats(const T *array) {
        for (int i = 0; i < size; i++) {
            feat[i] = array[i];
        }
    }

    template<typename T>
    void Vector<T>::print() {
        for (int i = 0; i < size; i++) {
            std::cout << feat[i] << " ";
        }
        std::cout << std::endl;
    }
}


#endif //SEGM_VECTOR_H
