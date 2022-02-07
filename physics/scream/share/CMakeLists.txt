
add_library(scream_share INTERFACE)
target_include_directories(scream_share INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(scream_share INTERFACE ekat)

