# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

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

cmake_minimum_required(VERSION 3.25)

# set the project name
project(tortoize VERSION 2.0.16 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(FindOrFetch)
include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(GenerateExportHeader)
include(CTest)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

# Optionally build a version to be installed inside CCP4
option(BUILD_FOR_CCP4 "Build a version to be installed in CCP4" OFF)

if(BUILD_FOR_CCP4)
	if("$ENV{CCP4}" STREQUAL "" OR NOT EXISTS $ENV{CCP4})
		message(FATAL_ERROR "tortoize: A CCP4 built was requested but CCP4 was not sourced")
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
	add_link_options(/NODEFAULTLIB:library)

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
find_package(mrc QUIET)

if(MRC_FOUND)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(STATUS "tortoize: Not using resources since mrc was not found")
	set(USE_RSRC OFF)
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

find_package(nlohmann_json REQUIRED)

find_or_fetch_package(mcfp VERSION 1.3.5 GIT_REPOSITORY https://github.com/mhekkel/libmcfp GIT_TAG b6c62a3)
find_or_fetch_package(cifpp VERSION 9 GIT_REPOSITORY https://github.com/PDB-REDO/libcifpp.git GIT_TAG v9.0.5 VARIABLES "CIFPP_SHARE_DIR")

if (NOT TARGET dssp)
	add_subdirectory(dssp EXCLUDE_FROM_ALL)
endif()

add_executable(tortoize
	${PROJECT_SOURCE_DIR}/src/tortoize.cpp
	${PROJECT_SOURCE_DIR}/src/tortoize-main.cpp)

if(USE_RSRC)
	list(APPEND RESOURCES
		${PROJECT_SOURCE_DIR}/rsrc/rama-data.bin
		${PROJECT_SOURCE_DIR}/rsrc/torsion-data.bin
		${CIFPP_SHARE_DIR}/mmcif_pdbx.dic
		${CIFPP_SHARE_DIR}/mmcif_ddl.dic
		${CIFPP_SHARE_DIR}/mmcif_ma.dic)

	mrc_target_resources(tortoize ${RESOURCES})
endif()

target_include_directories(tortoize PRIVATE ${PROJECT_BINARY_DIR})

target_link_libraries(tortoize dssp::dssp mcfp::mcfp cifpp::cifpp)

install(TARGETS tortoize
	RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)

if(NOT USE_RSRC)
	install(FILES ${PROJECT_SOURCE_DIR}/rsrc/rama-data.bin ${PROJECT_SOURCE_DIR}/rsrc/torsion-data.bin
		DESTINATION ${CIFPP_SHARE_DIR})
endif()

if(BUILD_WEBSERVICE)
	find_package(zeep 6 REQUIRED)

	if(NOT USE_RSRC)
		message(FATAL_ERROR "tortoize: Web service needs mrc/resources")
	endif()

	add_executable(tortoized
		${PROJECT_SOURCE_DIR}/src/tortoize.cpp
		${PROJECT_SOURCE_DIR}/src/tortoize-server-main.cpp
	)

	target_link_libraries(tortoized dssp::dssp mcfp::mcfp zeep::zeep)

	mrc_target_resources(tortoized ${RESOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/docroot/)

	install(TARGETS tortoized
		RUNTIME DESTINATION ${SBIN_INSTALL_DIR}
	)

	set(TORTOIZED_EXE "${CMAKE_INSTALL_FULL_SBINDIR}/tortoized")
	set(TORTOIZED_PID "/var/run/tortoized")

	configure_file(tortoized.service.in ${CMAKE_CURRENT_BINARY_DIR}/tortoized.service @ONLY)

	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tortoized.service" DESTINATION "/etc/systemd/system")
endif()

# manual
install(FILES doc/tortoize.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# Test code
if(BUILD_TESTING)
	add_subdirectory(test)
endif()
