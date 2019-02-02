#include "gaussiantest.h"
#include "imgproc/filter.h"
#include "test.h"

#include <cmath>

void gaussian_test()
{
    float expected[100] = {0.491836, 0.663096, 0.70131,  0.70131,  0.70131,  0.70131,  0.70131,  0.70131,  0.663096, 0.491836,
                           0.663096, 0.893991, 0.945511, 0.945511, 0.945511, 0.945511, 0.945511, 0.945511, 0.893991, 0.663096,
                           0.70131,  0.945511, 1,        1,        1,        1,        1,        1,        0.945511, 0.70131,
                           0.70131,  0.945511, 1,        1,        1,        1,        1,        1,        0.945511, 0.70131,
                           0.70131,  0.945511, 1,        1,        1,        1,        1,        1,        0.945511, 0.70131,
                           0.70131,  0.945511, 1,        1,        1,        1,        1,        1,        0.945511, 0.70131,
                           0.70131,  0.945511, 1,        1,        1,        1,        1,        1,        0.945511, 0.70131,
                           0.70131,  0.945511, 1,        1,        1,        1,        1,        1,        0.945511, 0.70131,
                           0.663096, 0.893991, 0.945511, 0.945511, 0.945511, 0.945511, 0.945511, 0.945511, 0.893991, 0.663096,
                           0.491836, 0.663096, 0.70131,  0.70131,  0.70131,  0.70131,  0.70131,  0.70131,  0.663096, 0.491836 };

    segm::Filter blur = segm::Filter::gaussian(5, 1.0);
    segm::Image<float> img(10, 10, 3);

    for (int i = 0; i < img.getHeight() * img.getWidth() * img.getBands(); i++)
        img(i) = 1.0f;

    segm::Image<float> res = segm::Filter::convolve(img, blur);

    for (int i = 0; i < res.getWidth(); i++) {
        for (int j = 0; j < res.getHeight(); j++) {
            for (int b = 0; b < res.getBands(); b++) {
                ASSERT_THROW((fabsf(res(i, j, b) - expected[i * 10 + j]) < 1e-4));
            }
        }
    }
}
