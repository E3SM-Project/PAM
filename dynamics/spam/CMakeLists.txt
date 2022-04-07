
add_library(dycore INTERFACE)
target_include_directories(dycore INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(src)
include_directories(src/driver)
include_directories(src/grids)
include_directories(src/models)
include_directories(src/operators)
include_directories(src/timesteppers)
include_directories(src/hamiltonians)

if ("${PAMC_MODEL}" STREQUAL "advection")
  include_directories(src/models/advection2D)
endif()

if ("${PAMC_MODEL}" STREQUAL "layer")
  include_directories(src/models/layermodel2D)
endif()

if ("${PAMC_MODEL}" STREQUAL "extruded")
  include_directories(src/models/extrudedmodel2D)
endif()
