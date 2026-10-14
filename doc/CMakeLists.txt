# manual 
find_program(PANDOC pandoc)
find_program(WEASYPRINT weasyprint)

if(PANDOC AND WEASYPRINT)
	file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pandoc-md.sh
		CONTENT "${PANDOC} ${CMAKE_CURRENT_SOURCE_DIR}/tortoize.1 -t markdown | sed -e \"s/\\\\\\' \\\\\\' (space)/  ' ' (space)/\" > ${CMAKE_CURRENT_SOURCE_DIR}/tortoize.md"
		FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)

	add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/tortoize.md
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tortoize.1
		COMMAND "${CMAKE_CURRENT_BINARY_DIR}/pandoc-md.sh"
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
		VERBATIM)

	add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/tortoize.html
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tortoize.1
		COMMAND ${PANDOC} -o tortoize.html tortoize.1
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
		VERBATIM)

	add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/tortoize.pdf
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tortoize.1
		COMMAND ${PANDOC} -o tortoize.pdf tortoize.1 -t html
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
		VERBATIM)

	add_custom_target(doc_files DEPENDS
		${CMAKE_CURRENT_SOURCE_DIR}/tortoize.html
		${CMAKE_CURRENT_SOURCE_DIR}/tortoize.pdf
		${CMAKE_CURRENT_SOURCE_DIR}/tortoize.md)

	add_dependencies(tortoize doc_files)
else()
	message(STATUS "Could not find pandoc, cannot recreate documentation files")
endif()
