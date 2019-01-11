#include "segmentation/watercuttest.h"

void watercut_test()
{
    float array[25] = {1, 1, 6, 5, 5,
                       1, 1, 5, 5, 5,
                       0, 1, 6, 6, 6,
                       0, 1, 9, 9, 9,
                       1, 1, 9, 9, 9};

    int m_array[25] = {-1, -1, -1, -1,  0,
                       -1, -1, -1, -1, -1,
                       -1, -1, -1, -1, -1,
                        1, -1, -1, -1, -1,
                       -1, -1, -1, -1,  2};

    segm::Image<int> markers(5, 5, m_array);
    segm::WaterCut cut(5, 5, 1, array);

    cut.run(markers);
    segm::Image<int> label = cut.getLabel();

    int true_label[25] = {1, 1, 0, 0, 0,
                          1, 1, 0, 0, 0,
                          1, 1, 0, 0, 0,
                          1, 1, 2, 2, 2,
                          1, 1, 2, 2, 2};

    for (int p = 0; p < 25; p++) {
        ASSERT_EQUAL(label(p), true_label[p]);
    }
}