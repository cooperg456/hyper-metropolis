##	---------------------------------------------------
##	--	Build Hyper Metropolis library: libmetro
##	---------------------------------------------------
add_library(libmetro)
add_library(hyper::metropolis ALIAS libmetro)

target_sources(libmetro 
PRIVATE 
	hamiltonian.cpp 
)

target_include_directories(libmetro
PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(libmetro 
PUBLIC 
	cxx_std_17
)

target_compile_definitions(libmetro 
PRIVATE 
	$<$<CONFIG:Debug>:DEBUG>
)

set_target_properties(libmetro 
PROPERTIES
	VERSION ${PROJECT_VERSION}
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES native
)
