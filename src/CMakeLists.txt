find_package(CBLAS REQUIRED)

set(SOURCE_FILES
        datatypes/forestingtransform.cpp
        datatypes/forestingtransform.h
        datatypes/heap.cpp
        datatypes/heap.h
        datatypes/image.cpp
        datatypes/image.h
        datatypes/matrix.cpp
        datatypes/matrix.h
        metriclearn/largemargin.cpp
        metriclearn/largemargin.h
        segmentation/watercut.cpp
        segmentation/watercut.h)

add_library(segm ${SOURCE_FILES})

target_include_directories(segm PRIVATE ${CBLAS_INCLUDES})
target_link_libraries(segm ${CBLAS_LIBRARIES})

target_include_directories(segm PUBLIC ${CMAKE_CURRENT_LIST_DIR})
