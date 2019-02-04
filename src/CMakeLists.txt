set(SOURCE_FILES
        datatypes/forestingtransform.cpp
        datatypes/forestingtransform.h
        datatypes/heap.cpp
        datatypes/heap.h
        datatypes/image.cpp
        datatypes/image.h
        datatypes/matrix.cpp
        datatypes/matrix.h
        datatypes/vector.cpp
        datatypes/vector.h
        imgproc/filter.cpp
        imgproc/filter.h
        math/linearalgebra.cpp
        math/linearalgebra.h
        math/statistics.cpp
        math/statistics.h
        metriclearn/largemargin.cpp
        metriclearn/largemargin.h
        metriclearn/neighborcompanalysis.cpp
        metriclearn/neighborcompanalysis.h
        segmentation/watercut.cpp
        segmentation/watercut.h
        utils/color.cpp
        utils/color.h)

add_library(segm ${SOURCE_FILES})

find_package(CBLAS REQUIRED)
target_include_directories(segm PRIVATE ${CBLAS_INCLUDES})
target_link_libraries(segm PRIVATE ${CBLAS_LIBRARIES})

find_package(LAPACK REQUIRED)
target_include_directories(segm PRIVATE ${LAPACK_INCLUDES})
target_link_libraries(segm PRIVATE ${LAPACK_LIBRARIES})

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