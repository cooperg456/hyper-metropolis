##	---------------------------------------------------
##	--	Build Hyper Metropolis test executables
##	---------------------------------------------------
set(TEST_SOURCES
    hamming.cpp
)

foreach(SRC_FILE IN LISTS TEST_SOURCES)
	get_filename_component(TARGET_NAME "${SRC_FILE}" NAME_WE)

	add_executable(${TARGET_NAME} ${SRC_FILE})

	target_link_libraries(${TARGET_NAME} PRIVATE hyper::metropolis)

	target_compile_features(${TARGET_NAME} PRIVATE cxx_std_17)

	target_compile_definitions(${TARGET_NAME} PRIVATE $<$<CONFIG:Debug>:DEBUG>)
endforeach()
