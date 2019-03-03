#include "datatypes/heaptest.h"
#include "imgproc/gaussiantest.h"
#include "math/pcatest.h"
#include "math/svdtest.h"
#include "metriclearn/largemargintest.h"
#include "metriclearn/ncatest.h"
#include "segmentation/dyntreetest.h"
#include "segmentation/watercuttest.h"
#include "utils/colortest.h"

int main()
{
    color_test();
    dyntree_test();
    gaussian_test();
    heap_test();
    largemargin_test();
    nca_test();
    pca_test();
    watercut_test();
}