#include "datatypes/heaptest.h"
#include "imgproc/gaussiantest.h"
#include "math/pcatest.h"
#include "math/svdtest.h"
#include "metriclearn/largemargintest.h"
#include "segmentation/watercuttest.h"
#include "utils/colortest.h"

int main()
{
    heap_test();
    gaussian_test();
    pca_test();
    svd_test();
    largemargin_test();
    watercut_test();
    color_test();
}