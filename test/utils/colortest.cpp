
#include "test.h"
#include "colortest.h"
#include "datatypes/image.h"

#include <cmath>

#include <iostream>

void color_test()
{
    int color[3] = {38, 50, 150};
    segm::Image<int> input(1, 1, 3, color);

    typedef segm::Image<double> Img;

    Img XYZ = input.convert(input.rgb, input.xyz, 255.0).convert(Img::xyz, Img::rgb);
    Img LAB = input.convert(input.rgb, input.lab, 255.0).convert(Img::lab, Img::rgb);
    Img YPbPr = input.convert(input.rgb, input.ypbpr, 255.0).convert(Img::ypbpr, Img::rgb);
    Img HSV = input.convert(input.rgb, input.hsv, 255.0).convert(Img::hsv, Img::rgb);

    for (int i = 0; i < input.getWidth(); i++) {
        for (int j = 0; j < input.getHeight(); j++) {
            for (int b = 0; b < input.getBands(); b++) {
                ASSERT_EQUAL((int) round(XYZ(i, j, b) * 255), input(i, j, b))
                ASSERT_EQUAL((int) round(LAB(i, j, b) * 255), input(i, j, b))
                ASSERT_EQUAL((int) round(YPbPr(i, j, b) * 255), input(i, j, b))
                ASSERT_EQUAL((int) round(HSV(i, j, b) * 255), input(i, j, b))
            }
        }
    }
}