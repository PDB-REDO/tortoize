CPMFindPackage(NAME Catch2 GIT_REPOSITORY "https://github.com/catchorg/Catch2" GIT_TAG v3.8.0 EXCLUDE_FROM_ALL YES)

add_executable(tortoize-unit-test
	${PROJECT_SOURCE_DIR}/test/tortoize-unit-test.cpp
	${PROJECT_SOURCE_DIR}/src/tortoize.cpp)

target_compile_definitions(tortoize-unit-test PUBLIC NOMINMAX=1)
target_include_directories(tortoize-unit-test PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR})

target_link_libraries(tortoize-unit-test nlohmann_json::nlohmann_json dssp::dssp Catch2::Catch2)

add_test(NAME tortoize-unit-test COMMAND $<TARGET_FILE:tortoize-unit-test> --data-dir ${PROJECT_SOURCE_DIR}/test)

if(USE_RSRC)
	mrc_target_resources(tortoize-unit-test ${RESOURCES})
endif()