
set(SCR_SHR_SRC pam_scream_routines.cpp)

add_library(scream_share ${SCR_SHR_SRC})
target_include_directories(scream_share PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(scream_share ekat)

