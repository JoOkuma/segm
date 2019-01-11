//
// Created by jookuma on 11/01/19.
//

#include "largemargin.h"

using namespace segm;

LargeMargin::LargeMargin(int _size, int dimension) :
    size(_size), d(dimension),
    L(d, d)
{
    for (int i = 0; i < d; i++) {
        L(i, i) = 1.0;
    }
}