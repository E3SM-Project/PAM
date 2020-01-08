
set(CPP_SRC 
  YAKL.cpp
  abcoefs.cpp
)

add_library(cpp_code STATIC ${CPP_SRC})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

