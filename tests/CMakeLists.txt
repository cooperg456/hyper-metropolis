##	---------------------------------------------------
##	--	Build Hyper Metropolis test executables
##	---------------------------------------------------
file(GLOB CONFIGURE_DEPENDS TEST_SOURCES "*.cpp")

foreach(SRC_FILE IN LISTS TEST_SOURCES)
	get_filename_component(TARGET_NAME "${SRC_FILE}" NAME_WE)

	add_executable(${TARGET_NAME} ${SRC_FILE})

	target_link_libraries(${TARGET_NAME} PRIVATE hyper::metropolis)

	target_compile_definitions(${TARGET_NAME} PRIVATE $<$<CONFIG:Debug>:DEBUG>)
endforeach()
