
set(SRC
 params.F90
 abcoefs.F90
)

add_library(fortran_code STATIC ${SRC})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

