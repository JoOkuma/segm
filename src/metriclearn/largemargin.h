//
// Created by jookuma on 11/01/19.
//

#ifndef SEGM_LARGEMARGIN_H
#define SEGM_LARGEMARGIN_H

#include "datatypes/matrix.h"

namespace segm
{
    class LargeMargin
    {
    public:
        LargeMargin(int _size, int dimension);
        ~LargeMargin();

    private:
        Matrix<float> L;

        int size;
        int d;
    };
}


#endif //SEGM_LARGEMARGIN_H
