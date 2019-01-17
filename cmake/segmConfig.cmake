add_library(segm STATIC IMPORTED)
find_library(SEGM_LIBRARY_PATH segm HINTS "${CMAKE_CURRENT_LIST_DIR}/../../")
set_target_properties(segm PROPERTIES IMPORTED_LOCATION "${SEGM_LIBRARY_PATH}")