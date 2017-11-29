################################################################################
# Project:  NextGIS command line utilities projects
# Purpose:  CMake build scripts
# Author:   Dmitry Baryshnikov, dmitry.baryshnikov@nexgis.com
################################################################################
# Copyright (C) 2017, NextGIS <info@nextgis.com>
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this script.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

cmake_minimum_required (VERSION 2.8.12 FATAL_ERROR)

project(ng_utilities)

# Some init settings
set(CMAKE_COLOR_MAKEFILE ON)
# Set path to additional CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

if(OSX_FRAMEWORK AND (BUILD_SHARED_LIBS OR BUILD_STATIC_LIBS))
  message(FATAL_ERROR "Only OSX_FRAMEWORK key or any or both BUILD_SHARED_LIBS
                       and BUILD_STATIC_LIBS keys are permitted")
endif()

if(OSX_FRAMEWORK)
  set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE INTERNAL "Installation directory for executables" FORCE)
  set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/Library/Frameworks" CACHE INTERNAL "Installation directory for libraries" FORCE)
  set(SKIP_INSTALL_EXECUTABLES OFF)
  set(SKIP_INSTALL_EXPORT ON)
  set(SKIP_INSTALL_DOCS ON)
  set(CMAKE_MACOSX_RPATH ON)
  set(GDAL_BUILD_DOCS OFF CACHE INTERNAL "Don't build docs on Mac OS X" FORCE)
else()
    include(GNUInstallDirs)

    set(INSTALL_BIN_DIR ${CMAKE_INSTALL_BINDIR} CACHE INTERNAL "Installation directory for executables" FORCE)
    set(INSTALL_LIB_DIR ${CMAKE_INSTALL_LIBDIR} CACHE INTERNAL "Installation directory for libraries" FORCE)
    set(INSTALL_INC_DIR "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}" CACHE INTERNAL "Installation directory for headers" FORCE)
    set(INSTALL_SHARE_DIR "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}" CACHE INTERNAL "Installation directory for share data" FORCE)
    set(INSTALL_PKGCONFIG_DIR "${INSTALL_LIB_DIR}/pkgconfig" CACHE INTERNAL "Installation directory for pkgconfig (.pc) files" FORCE)
    set(INSTALL_MAN_DIR ${CMAKE_INSTALL_MANDIR} CACHE INTERNAL "Installation directory for manual pages" FORCE)
    set(INSTALL_DOC_DIR ${CMAKE_INSTALL_DOCDIR} CACHE INTERNAL "Installation directory for doc pages" FORCE)
endif()

include(util)

check_version(NG_MAJOR_VERSION NG_MINOR_VERSION NG_REV_VERSION)
set(VERSION "${NG_MAJOR_VERSION}.${NG_MINOR_VERSION}.${NG_REV_VERSION}")
report_version(${PROJECT_NAME} ${VERSION})

include(FindAnyProject)
find_anyproject(GDAL REQUIRED)
find_anyproject(GEOS REQUIRED)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

set(INSTALL_TARGETS)

macro(ADD_EXE target sources)
    add_executable(${target} ${sources})
    target_link_libraries(${target} ${TARGET_LINK_LIB})
    if(OSX_FRAMEWORK)
        set_target_properties(${target} PROPERTIES INSTALL_RPATH "@executable_path/../../Library/Frameworks")
    endif()
    # set_target_properties(${target} PROPERTIES CXX_STANDARD 11)
    set(INSTALL_TARGETS ${INSTALL_TARGETS} ${target})
endmacro()

set(NG_CUTTER_H
    commonutils.h
    ogr2ogr_lib.h
    ogrlayerdecorator.h
)

set(NG_CUTTER_SRC
    ng_cutter.cpp
    commonutils.cpp
    ogr2ogr_lib.cpp
    ogrlayerdecorator.cpp
)

set(SRC ${NG_CUTTER_SRC} ${NG_CUTTER_H})
add_exe(ng_cutter "${SRC}")

if(NOT SKIP_INSTALL_EXECUTABLES AND NOT SKIP_INSTALL_ALL )
    install(TARGETS ${INSTALL_TARGETS}
    #    EXPORT ${INSTALL_TARGETS}
        RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT applications
        ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT applications
        LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT applications
        FRAMEWORK DESTINATION ${INSTALL_BIN_DIR}
    )
endif()
