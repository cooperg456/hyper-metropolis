##	---------------------------------------------------
##	--	Build Hyper Metropolis library: libmetro
##	---------------------------------------------------
file(GLOB CPP_FILES "${CMAKE_CURRENT_LIST_DIR}/*.cpp")
file(GLOB HPP_FILES "${CMAKE_CURRENT_LIST_DIR}/*.hpp")

add_library(metro)

target_sources(metro PRIVATE ${CPP_FILES} ${HPP_FILES} PUBLIC "${CMAKE_SOURCE_DIR}/include/hyper-metropolis/metro.hpp")

target_include_directories(metro PRIVATE "${CMAKE_SOURCE_DIR}/dependencies/vma/include")

find_package(Vulkan REQUIRED)

target_link_libraries(metro PRIVATE Vulkan::Vulkan)

target_compile_definitions(metro PRIVATE $<$<CONFIG:Debug>:DEBUG>)

set_target_properties(metro PROPERTIES VERSION ${PROJECT_VERSION})