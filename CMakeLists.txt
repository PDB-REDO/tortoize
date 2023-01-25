# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.15)

# set the project name
project(tortoize VERSION 2.0.11 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(Dart)
include(FindFilesystem)
include(GenerateExportHeader)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Filesystem REQUIRED)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

# Optionally build a version to be installed inside CCP4
option(BUILD_FOR_CCP4 "Build a version to be installed in CCP4" OFF)

if(BUILD_FOR_CCP4)
	set(USE_RSRC 0)
	if("$ENV{CCP4}" STREQUAL "" OR NOT EXISTS $ENV{CCP4})
		message(FATAL_ERROR "A CCP4 built was requested but CCP4 was not sourced")
	else()
		list(APPEND CMAKE_MODULE_PATH "$ENV{CCP4}")
		list(APPEND CMAKE_PREFIX_PATH "$ENV{CCP4}")
		set(CMAKE_INSTALL_PREFIX "$ENV{CCP4}")

		if(WIN32)
			set(BUILD_SHARED_LIBS ON)
		endif()
	endif()
endif()

if(MSVC)
	# make msvc standards compliant...
	add_compile_options(/permissive-)

	macro(get_WIN32_WINNT version)
		if(WIN32 AND CMAKE_SYSTEM_VERSION)
			set(ver ${CMAKE_SYSTEM_VERSION})
			string(REPLACE "." "" ver ${ver})
			string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})

			set(${version} "0x${ver}")
		endif()
	endmacro()

	get_WIN32_WINNT(ver)
	add_definitions(-D_WIN32_WINNT=${ver})
endif()

# Create a revision file, containing the current git version info, if any
include(VersionString)
write_version_header(${PROJECT_SOURCE_DIR}/src)

# Optionally use mrc to create resources
if(NOT DEFINED USE_RSRC)
	find_package(Mrc QUIET)

	if(MRC_FOUND)
		option(USE_RSRC "Use mrc to create resources" ON)
	else()
		message(STATUS "Not using resources since mrc was not found")
		set(USE_RSRC OFF)
	endif()
endif()

if(USE_RSRC)
	set(USE_RSRC 1)
	add_compile_definitions(USE_RSRC)
endif()

# Optionally build a webservice
option(BUILD_WEBSERVICE "Build a version with a webservice daemon" OFF)

# libraries
set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

if(NOT PDB_REDO_META)
	find_package(libmcfp REQUIRED)
	find_package(dssp REQUIRED)
	find_package(cifpp REQUIRED)
	find_package(zeep 5.1.8 REQUIRED)
endif()

add_executable(tortoize
	${PROJECT_SOURCE_DIR}/src/tortoize.cpp
	${PROJECT_SOURCE_DIR}/src/tortoize-main.cpp
	${TORTOIZE_RESOURCE})

if(BUILD_WEBSERVICE)
	target_compile_definitions(tortoize PRIVATE WEBSERVICE)
endif()

if(USE_RSRC)
	list(APPEND RESOURCES
		${PROJECT_SOURCE_DIR}/rsrc/rama-data.bin
		${PROJECT_SOURCE_DIR}/rsrc/torsion-data.bin
		${CIFPP_SHARE_DIR}/mmcif_pdbx.dic
		${CIFPP_SHARE_DIR}/mmcif_ddl.dic
		${CIFPP_SHARE_DIR}/mmcif_ma.dic)

	if(BUILD_WEBSERVICE)
		list(APPEND RESOURCES ${PROJECT_SOURCE_DIR}/docroot/)
	endif()

	mrc_target_resources(tortoize ${RESOURCES})
endif()

target_include_directories(tortoize PRIVATE ${PROJECT_BINARY_DIR})

if(PDB_REDO_META)
	target_include_directories(tortoize PRIVATE ${CMAKE_PROJECT_DIR}/dssp/include)
endif()

target_link_libraries(tortoize dssp::dssp cifpp::cifpp zeep::zeep std::filesystem libmcfp::libmcfp)
target_compile_definitions(tortoize PUBLIC NOMINMAX=1)

install(TARGETS ${PROJECT_NAME}
	RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)

if(NOT USE_RSRC)
	install(FILES ${PROJECT_SOURCE_DIR}/rsrc/rama-data.bin ${PROJECT_SOURCE_DIR}/rsrc/torsion-data.bin
		DESTINATION ${CIFPP_SHARE_DIR})
endif()

# manual
install(FILES doc/tortoize.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# Test code
if(ENABLE_TESTING)
	add_executable(tortoize-unit-test
		${PROJECT_SOURCE_DIR}/test/tortoize-unit-test.cpp
		${PROJECT_SOURCE_DIR}/src/tortoize.cpp)

	target_compile_definitions(tortoize-unit-test PUBLIC NOMINMAX=1)
	target_include_directories(tortoize-unit-test PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR})

	target_link_libraries(tortoize-unit-test dssp::dssp cifpp::cifpp zeep::zeep std::filesystem)

	add_test(NAME tortoize-unit-test COMMAND $<TARGET_FILE:tortoize-unit-test> -- ${PROJECT_SOURCE_DIR}/test)

	if(USE_RSRC)
		mrc_target_resources(tortoize-unit-test ${RESOURCES})
	endif()
endif()
