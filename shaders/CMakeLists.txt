##	---------------------------------------------------
##	--	Build HLSL shaders to SPIR-V
##	---------------------------------------------------
file(GLOB HLSL_FILES "${CMAKE_CURRENT_LIST_DIR}/*.hlsl")

set(SPV_FILES "")

foreach(HLSL_FILE ${HLSL_FILES}) 
	get_filename_component(HLSL_FILENAME ${HLSL_FILE} NAME_WLE) 
	set(OUTPUT_FILE "${CMAKE_BINARY_DIR}/shaders/${HLSL_FILENAME}.spv")
	add_custom_command( 
		OUTPUT ${OUTPUT_FILE}
		COMMAND glslangValidator 
				-S comp 
				$<$<CONFIG:Debug>:-g> 
				$<$<CONFIG:Debug>:-H>
				-o ${OUTPUT_FILE} 
				-V 
				-D ${HLSL_FILE}
		DEPENDS ${HLSL_FILE}
		COMMENT "Building shader ${OUTPUT_FILE}"
	) 
	list(APPEND SPV_FILES ${OUTPUT_FILE})
endforeach()

add_custom_target(shaders ALL DEPENDS ${SPV_FILES})