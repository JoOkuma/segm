
#include "test.h"
#include "segmentation/isf.h"
#include "segmentation/isftest.h"

void isf_test()
{
    float array[32] = {  1,   1,   1,   1,  64,  64,  64,  64,
                         1,   1,   1, 255,  64,  64,  64, 127,
                         1, 255, 255, 255, 127, 127,  64, 127,
                       255, 255, 255, 255, 127, 127, 127, 127 };

    segm::Image<float> original(8, 4, 1, array);
    segm::ISF forest(original);

    forest.run(4, 2);
    segm::Image<int> label = forest.getLabel();

    int true_label[32] = { 0, 0, 0, 0, 2, 2, 2, 2,
                           0, 0, 0, 1, 2, 2, 2, 3,
                           0, 1, 1, 1, 3, 3, 2, 3,
                           1, 1, 1, 1, 3, 3, 3, 3 };

    for (int p = 0; p < label.getSize(); p++) {
        ASSERT_EQUAL(label(p), true_label[p]);
    }
}