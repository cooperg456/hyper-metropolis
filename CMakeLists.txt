cmake_minimum_required(VERSION 3.23)

project(HYPER_METROPOLIS VERSION 0.0.1 LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

option(METRO_BUILD_TESTS "Build test programs" ON)

# ---------------------------------------------------
# -- header only library
# ---------------------------------------------------

add_library(metrolib INTERFACE)
add_library(hyper::metropolis ALIAS metrolib)

target_include_directories(metrolib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_features(metrolib INTERFACE cxx_std_17)

##	---------------------------------------------------
##	--	add subdirectories
##	---------------------------------------------------

if (METRO_BUILD_TESTS)
	add_subdirectory(tests)
endif()
