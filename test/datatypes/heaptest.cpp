#include "test.h"
#include "datatypes/heap.h"
#include "datatypes/heaptest.h"

void heap_test()
{
    float input[6] = {10, 20, 30, 5, 15, 99};
    segm::Heap heap(2, 3);
    heap.setValues(input);

    for (int i = 0; i < 6; i++)
        heap.insert(i);

    int min_heap[6] = {3, 0, 4, 1, 2, 5};
    for (int i = 0; i < 6; i++) {
        int p = heap.pop();
        ASSERT_EQUAL(p, min_heap[i]);
    }

    heap.setPolicy(segm::Heap::max_value);
    heap.reset();

    for (int i = 0; i < 6; i++)
        heap.insert(i);

    int max_heap[6] = {5, 2, 1, 4, 0, 3};
    for (int i = 0; i < 6; i++) {
        int p = heap.pop();
        ASSERT_EQUAL(p, max_heap[i]);
    }
}