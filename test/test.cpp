#include "datatypes/heaptest.h"
#include "imgproc/gaussiantest.h"
#include "math/pcatest.h"
#include "math/svdtest.h"
#include "metriclearn/largemargintest.h"
#include "metriclearn/ncatest.h"
#include "segmentation/watercuttest.h"
#include "utils/colortest.h"

int main()
{
    heap_test();
    gaussian_test();
    pca_test();
    largemargin_test();
    nca_test();
    watercut_test();
    color_test();
}