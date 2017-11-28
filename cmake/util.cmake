################################################################################
# Project:  external projects
# Purpose:  CMake build scripts
# Author:   Dmitry Baryshnikov, polimax@mail.ru
################################################################################
# Copyright (C) 2015, NextGIS <info@nextgis.com>
# Copyright (C) 2015 Dmitry Baryshnikov
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


function(check_version major minor rev)

    # parse the version number from gdal_version.h and include in
    # major, minor and rev parameters
    set(VERSION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/version.h)

    file(READ ${VERSION_FILE} _VERSION_H_CONTENTS)

    string(REGEX MATCH "NGU_VERSION_MAJOR[ \t]+([0-9]+)"
      NGU_MAJOR_VERSION ${_VERSION_H_CONTENTS})
    string (REGEX MATCH "([0-9]+)"
      NGU_MAJOR_VERSION ${NGU_MAJOR_VERSION})
    string(REGEX MATCH "NGU_VERSION_MINOR[ \t]+([0-9]+)"
      NGU_MINOR_VERSION ${_VERSION_H_CONTENTS})
    string (REGEX MATCH "([0-9]+)"
      NGU_MINOR_VERSION ${NGU_MINOR_VERSION})
    string(REGEX MATCH "NGU_VERSION_REV[ \t]+([0-9]+)"
      NGU_REV_VERSION ${_VERSION_H_CONTENTS})
    string (REGEX MATCH "([0-9]+)"
      NGU_REV_VERSION ${NGU_REV_VERSION})

    set(${major} ${NGU_MAJOR_VERSION} PARENT_SCOPE)
    set(${minor} ${NGU_MINOR_VERSION} PARENT_SCOPE)
    set(${rev}   ${NGU_REV_VERSION}   PARENT_SCOPE)

    # Store version string in file for installer needs
    file(TIMESTAMP ${VERSION_FILE} VERSION_DATETIME "%Y-%m-%d %H:%M:%S" UTC)
    file(WRITE ${CMAKE_BINARY_DIR}/version.str "${NGU_MAJOR_VERSION}.${NGU_MINOR_VERSION}.${NGU_REV_VERSION}\n${VERSION_DATETIME}")

endfunction(check_version)

function(report_version name ver)

    string(ASCII 27 Esc)
    set(BoldYellow  "${Esc}[1;33m")
    set(ColourReset "${Esc}[m")

    message(STATUS "${BoldYellow}${name} version ${ver}${ColourReset}")

endfunction()


# macro to find packages on the host OS
macro( find_exthost_package )
    if(CMAKE_CROSSCOMPILING)
        set( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER )
        set( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY NEVER )
        set( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE NEVER )

        find_package( ${ARGN} )

        set( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM ONLY )
        set( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY )
        set( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY )
    else()
        find_package( ${ARGN} )
    endif()
endmacro()


# macro to find programs on the host OS
macro( find_exthost_program )
    if(CMAKE_CROSSCOMPILING)
        set( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER )
        set( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY NEVER )
        set( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE NEVER )

        find_program( ${ARGN} )

        set( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM ONLY )
        set( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY )
        set( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY )
    else()
        find_program( ${ARGN} )
    endif()
endmacro()


# macro to find path on the host OS
macro( find_exthost_path )
    if(CMAKE_CROSSCOMPILING)
        set( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER )
        set( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY NEVER )
        set( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE NEVER )

        find_path( ${ARGN} )

        set( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM ONLY )
        set( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY )
        set( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY )
    else()
        find_path( ${ARGN} )
    endif()
endmacro()
