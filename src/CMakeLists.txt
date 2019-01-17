set(SOURCE_FILES
        datatypes/forestingtransform.cpp
        datatypes/heap.cpp
        datatypes/image.cpp
        datatypes/matrix.cpp
        metriclearn/largemargin.cpp
        segmentation/watercut.cpp)

set(HEADER_FILES
        datatypes/forestingtransform.h
        datatypes/heap.h
        datatypes/image.h
        datatypes/matrix.h
        metriclearn/largemargin.h
        segmentation/watercut.h)

add_library(segm ${SOURCE_FILES} ${HEADER_FILES})

find_package(CBLAS REQUIRED)
target_include_directories(segm PRIVATE ${CBLAS_INCLUDES})
target_link_libraries(segm PRIVATE ${CBLAS_LIBRARIES})

find_package(OpenMP)
target_link_libraries(segm PRIVATE OpenMP::OpenMP_CXX)

target_include_directories(segm PUBLIC ${CMAKE_CURRENT_LIST_DIR})

install(FILES ${HEADER_FILES} DESTINATION include)
install(TARGETS segm DESTINATION lib)
install(FILES ../cmake/segmConfig.cmake
        DESTINATION lib/cmake/segm)