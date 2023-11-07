project(map C CXX)
set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "0")
cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR) 

add_library(${PROJECT_NAME} STATIC "")

set(CMAKE_C_STANDARD 99)

target_sources(
	${PROJECT_NAME}
	PRIVATE
		"${CMAKE_CURRENT_LIST_DIR}/map.c"
)

target_include_directories(
	${PROJECT_NAME}
	PUBLIC
		"${CMAKE_CURRENT_LIST_DIR}"
)
