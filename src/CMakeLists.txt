set(SOURCE_FILES
        datatypes/image.cpp
        datatypes/image.h
        datatypes/heap.cpp
        datatypes/heap.h
        datatypes/matrix.cpp
        datatypes/matrix.h
        metriclearn/largemargin.cpp
        metriclearn/largemargin.h
        segmentation/watercut.cpp
        segmentation/watercut.h)

add_library(segm ${SOURCE_FILES})
target_include_directories(segm PUBLIC ${CMAKE_CURRENT_LIST_DIR})
