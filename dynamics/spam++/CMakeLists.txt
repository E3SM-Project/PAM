


cmake_minimum_required(VERSION 3.0)
project(spam++ CXX)

if ("${GENERAL_COMPILE_CONST}" STREQUAL "")
file(READ src/compile_consts.build general_defines)
else()
file(READ ${GENERAL_COMPILE_CONST} general_defines)
endif()
# Convert file contents into a CMake list (where each element in the list
# is one line of the file)
#
STRING(REGEX REPLACE ";" "\\\\;" general_defines "${general_defines}")
STRING(REGEX REPLACE "\n" ";" general_defines "${general_defines}")

foreach(definecmd ${general_defines})
if ("${definecmd}" MATCHES "\#")
else ()
message("${definecmd}")
add_definitions(${definecmd})
endif()
endforeach()







# basically the build script takes 
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS}")

include_directories(${PNETCDF_INCLUDE})
link_directories(${PNETCDF_LIB})
link_libraries(pnetcdf)

add_subdirectory(externals/YAKL)
include_directories(externals/YAKL)
include_directories(externals/cub)
include_directories(externals/hipCUB)
include_directories(externals/rocPRIM)

include_directories(src)
include_directories(src/driver)
include_directories(src/grids)
include_directories(src/models)
include_directories(src/operators)
include_directories(src/timesteppers)
include_directories(src/hamiltonians)

if ("${MODEL}" STREQUAL "advection2D")
  include_directories(src/models/advection2D)
  add_executable(advection2D driver/driver.cpp)
  set_target_properties(advection2D PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(advection2D yakl pnetcdf)
  set_source_files_properties(driver.cpp PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${YAKL_CXX_FLAGS}")

  if ("${MODEL_COMPILE_CONST}" STREQUAL "")
  file(READ src/models/advection2D/model_consts.build model_defines)
  else()
  file(READ ${MODEL_COMPILE_CONST} model_defines)
  endif()
  # Convert file contents into a CMake list (where each element in the list
  # is one line of the file)
  #
  STRING(REGEX REPLACE ";" "\\\\;" model_defines "${model_defines}")
  STRING(REGEX REPLACE "\n" ";" model_defines "${model_defines}")

  foreach(definecmd ${model_defines})
  if ("${definecmd}" MATCHES "\#")
  else ()
  message("${definecmd}")
  add_definitions(${definecmd})
  endif()
  endforeach()

endif()

if ("${MODEL}" STREQUAL "extrudedadvection2D")

  include_directories(src/models/extrudedadvection2D)
  add_executable(extrudedadvection2D driver/driver.cpp)
  set_target_properties(extrudedadvection2D PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(extrudedadvection2D yakl pnetcdf)
  set_source_files_properties(driver.cpp PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${YAKL_CXX_FLAGS}")

  if ("${MODEL_COMPILE_CONST}" STREQUAL "")
  file(READ src/models/extrudedadvection2D/model_consts.build model_defines)
  else()
  file(READ ${MODEL_COMPILE_CONST} model_defines)
  endif()
  # Convert file contents into a CMake list (where each element in the list
  # is one line of the file)
  #
  STRING(REGEX REPLACE ";" "\\\\;" model_defines "${model_defines}")
  STRING(REGEX REPLACE "\n" ";" model_defines "${model_defines}")

  foreach(definecmd ${model_defines})
  if ("${definecmd}" MATCHES "\#")
  else ()
  message("${definecmd}")
  add_definitions(${definecmd})
  endif()
  endforeach()

  add_definitions(-D_EXTRUDED)

endif()


if ("${MODEL}" STREQUAL "tswe2D" OR "${MODEL}" STREQUAL "swe2D" OR "${MODEL}" STREQUAL "ce2D" OR "${MODEL}" STREQUAL "mce2D")
  include_directories(src/models/layermodel2D)
  add_executable(layermodel2D driver/driver.cpp)
  set_target_properties(layermodel2D PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(layermodel2D yakl pnetcdf)
  set_source_files_properties(driver.cpp PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${YAKL_CXX_FLAGS}")

  if ("${MODEL_COMPILE_CONST}" STREQUAL "")
  file(READ src/models/layermodel2D/model_consts.build model_defines)
  else()
  file(READ ${MODEL_COMPILE_CONST} model_defines)
  endif()
  # Convert file contents into a CMake list (where each element in the list
  # is one line of the file)
  #
  STRING(REGEX REPLACE ";" "\\\\;" model_defines "${model_defines}")
  STRING(REGEX REPLACE "\n" ";" model_defines "${model_defines}")

  foreach(definecmd ${model_defines})
  if ("${definecmd}" MATCHES "\#")
  else ()
  message("${definecmd}")
  add_definitions(${definecmd})
  endif()
  endforeach()

  if ("${MODEL}" STREQUAL "swe2D")
  add_definitions(-D_SWE)
  endif()

  if ("${MODEL}" STREQUAL "tswe2D")
  add_definitions(-D_TSWE)
  endif()

  if ("${MODEL}" STREQUAL "ce2D")
  add_definitions(-D_CE)
  endif()

  if ("${MODEL}" STREQUAL "mce2D")
  add_definitions(-D_MCE)
  endif()

endif()

if ("${MODEL}" STREQUAL "tswe2Dext" OR "${MODEL}" STREQUAL "swe2Dext" OR "${MODEL}" STREQUAL "ce2Dext" OR "${MODEL}" STREQUAL "mce2Dext")
  include_directories(src/models/extrudedmodel2D)
  add_executable(extrudedmodel2D driver/driver.cpp)
  set_target_properties(extrudedmodel2D PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(extrudedmodel2D yakl pnetcdf)
  set_source_files_properties(driver.cpp PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${YAKL_CXX_FLAGS}")

  if ("${MODEL_COMPILE_CONST}" STREQUAL "")
  file(READ src/models/extrudedmodel2D/model_consts.build model_defines)
  else()
  file(READ ${MODEL_COMPILE_CONST} model_defines)
  endif()
  # Convert file contents into a CMake list (where each element in the list
  # is one line of the file)
  #
  STRING(REGEX REPLACE ";" "\\\\;" model_defines "${model_defines}")
  STRING(REGEX REPLACE "\n" ";" model_defines "${model_defines}")

  foreach(definecmd ${model_defines})
  if ("${definecmd}" MATCHES "\#")
  else ()
  message("${definecmd}")
  add_definitions(${definecmd})
  endif()
  endforeach()

  add_definitions(-D_EXTRUDED)

  if ("${MODEL}" STREQUAL "swe2Dext")
  add_definitions(-D_SWE)
  endif()

  if ("${MODEL}" STREQUAL "tswe2Dext")
  add_definitions(-D_TSWE)
  endif()

  if ("${MODEL}" STREQUAL "ce2Dext")
  add_definitions(-D_CE)
  endif()

  if ("${MODEL}" STREQUAL "mce2Dext")
  add_definitions(-D_MCE)
  endif()

endif()
