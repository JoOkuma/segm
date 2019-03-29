add_executable(main_test
        test.cpp
        datatypes/heaptest.cpp
        datatypes/heaptest.h
        imgproc/gaussiantest.cpp
        imgproc/gaussiantest.h
        math/pcatest.cpp
        math/pcatest.h
        metriclearn/largemargintest.cpp
        metriclearn/largemargintest.h
        metriclearn/ncatest.cpp
        metriclearn/ncatest.h
        segmentation/dyntreetest.cpp
        segmentation/dyntreetest.h
        segmentation/isftest.cpp
        segmentation/isftest.h
        segmentation/watercuttest.cpp
        segmentation/watercuttest.h
        utils/colortest.cpp
        utils/colortest.h)

target_include_directories(main_test PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(main_test segm)

add_test(test main_test)