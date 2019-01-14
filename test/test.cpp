#include "datatypes/heaptest.h"
#include "metriclearn/largemargintest.h"
#include "segmentation/watercuttest.h"

int main(void)
{
    heap_test();
    largemargin_test();
    watercut_test();
}