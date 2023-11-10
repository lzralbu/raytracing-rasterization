cmake_minimum_required(VERSION 3.15)
project(computer-graphics)

set(CMAKE_C_STANDARD 99)

find_package(raylib 4.0 REQUIRED)

add_executable(${PROJECT_NAME} src/main.c)

add_subdirectory(libs)

target_compile_options(${PROJECT_NAME} PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>: /W4 /Zc:preprocessor>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -Wall -Wextra -Wpedantic>
)

target_link_libraries(${PROJECT_NAME} raylib)
target_link_libraries(${PROJECT_NAME} xml)
target_link_libraries(${PROJECT_NAME} map)
target_link_libraries(${PROJECT_NAME} vec)

# Checks if OSX and links appropriate frameworks (only required on MacOS)
if (APPLE)
    target_link_libraries(${PROJECT_NAME} "-framework IOKit")
    target_link_libraries(${PROJECT_NAME} "-framework Cocoa")
    target_link_libraries(${PROJECT_NAME} "-framework OpenGL")
endif()

# add_dependencies(${PROJECT_NAME} copy_assets)

# add_custom_target(copy_assets
#     COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/copy-assets.cmake
# )
