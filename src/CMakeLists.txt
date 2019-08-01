set(SOURCE_FILES
        datatypes/forestingtransform.cpp
        datatypes/forestingtransform.h
        datatypes/graph.cpp
        datatypes/graph.h
        datatypes/heap.cpp
        datatypes/heap.h
        datatypes/image.cpp
        datatypes/image.h
        datatypes/optimumpathforest.cpp
        datatypes/optimumpathforest.h
        imgproc/filter.cpp
        imgproc/filter.h
        math/statistics.cpp
        math/statistics.h
        metriclearn/largemargin.cpp
        metriclearn/largemargin.h
        metriclearn/neighborcompanalysis.cpp
        metriclearn/neighborcompanalysis.h
        segmentation/closestdyntree.cpp
        segmentation/closestdyntree.h
        segmentation/dtisf.cpp
        segmentation/dtisf.h
        segmentation/dyntree.cpp
        segmentation/dyntree.h
        segmentation/watercut.cpp
        segmentation/watercut.h
        segmentation/isf.cpp
        segmentation/isf.h
        utils/algorithm.cpp
        utils/algorithm.h
        utils/color.cpp
        utils/color.h)

add_library(segm ${SOURCE_FILES})

find_package(OpenMP)
if (OPENMP_FOUND)
    target_link_libraries(segm PRIVATE OpenMP::OpenMP_CXX)
endif()

target_include_directories(segm PUBLIC ${CMAKE_CURRENT_LIST_DIR})

install(TARGETS segm DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/
        DESTINATION include
        FILES_MATCHING PATTERN "*/*.h")
install(FILES ../cmake/segmConfig.cmake
        DESTINATION lib/cmake/segm)