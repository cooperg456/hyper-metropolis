cmake_minimum_required(VERSION 3.23)

project(HYPER_METROPOLIS VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

string(COMPARE EQUAL "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" METRO_STANDALONE)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(METRO_BUILD_EXAMPLES "Build example programs" ${METRO_STANDALONE})
option(METRO_BUILD_TESTS "Build test programs" ${METRO_STANDALONE})
option(METRO_BUILD_SHADERS "Build shader binaries" ON)
option(METRO_BUILD_PYTHON "Build python bindings" OFF)
option(METRO_BUILD_DOCS "Build documentations" OFF)

##	---------------------------------------------------
##	--	Add subdirectories
##	---------------------------------------------------
add_subdirectory(src)

if (METRO_BUILD_EXAMPLES)
	add_subdirectory(examples)
endif()

if (METRO_BUILD_TESTS)
	add_subdirectory(tests)
endif()

if (METRO_BUILD_SHADERS)
	add_subdirectory(shaders)
endif()

if (METRO_BUILD_PYTHON)
	add_subdirectory(python)
endif()

if (METRO_BUILD_DOCS)
	add_subdirectory(docs)
endif()