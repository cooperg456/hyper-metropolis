cmake_minimum_required(VERSION 3.23)

project(HYPER_METROPOLIS VERSION 0.0.1 LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

string(COMPARE EQUAL "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" METRO_STANDALONE)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(METRO_BUILD_TESTS "Build test programs" ${METRO_STANDALONE})

##	---------------------------------------------------
##	--	Add subdirectories
##	---------------------------------------------------
add_subdirectory(src)

if (METRO_BUILD_TESTS)
	add_subdirectory(tests)
endif()
