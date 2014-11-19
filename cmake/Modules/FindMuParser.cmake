# - Find MuParser
# Downloaded from
#  http://trac.openturns.org/export/3306/openturns/trunk/cmake/FindMuParser.cmake
# on November 10, 2014
#
# muParser is an extensible high performance math expression parser library written in C++
# http://muparser.sourceforge.net
#
# The module defines the following variables:
#  MUPARSER_FOUND        - True if MuParser found.
#  MUPARSER_INCLUDE_DIRS - where to find muParser.h, etc.
#  MUPARSER_LIBRARIES    - List of libraries when using MuParser.
#
#=============================================================================
# Note: the license on this file was originally identified as "the OSI-approved
# BSD license" without specifying which BSD license was intended. The following
# is the text of the BSD 3-clause license, which is the more restrictive of the
# two options.
#=============================================================================
# Copyright (C) 2005-2013 EDF-EADS-Phimeca
#
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this
# list of conditions and the following disclaimer in the documentation and/or other
# materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may
# be used to endorse or promote products derived from this software without specific
# prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

if ( MUPARSER_INCLUDE_DIR AND MUPARSER_LIBRARIES )
  # Already in cache, be silent
  set ( MuParser_FIND_QUIETLY TRUE )
endif ()

find_path ( MUPARSER_INCLUDE_DIR muParser.h
            PATHS
            C:/muparser/include
            PATH_SUFFIXES muParser )
            
set ( MUPARSER_INCLUDE_DIRS ${MUPARSER_INCLUDE_DIR} )

# version
set ( _VERSION_FILE ${MUPARSER_INCLUDE_DIR}/muParserDef.h )
if ( EXISTS ${_VERSION_FILE} )
  file ( STRINGS ${_VERSION_FILE} _MUPARSER_VERSION_LINE REGEX "#define[ ]+MUP_VERSION[ ]+_T\\(\".*\"\\)" )
  if ( _MUPARSER_VERSION_LINE )
    string ( REGEX REPLACE ".*_T\\(\"(.*)\"\\)" "\\1" MUPARSER_VERSION_STRING ${_MUPARSER_VERSION_LINE} )
    if ( MUPARSER_VERSION_STRING )
      string ( REGEX REPLACE "([0-9]+)\\..*" "\\1" MUPARSER_MAJOR_VERSION ${MUPARSER_VERSION_STRING} )
      string ( REGEX REPLACE "[0-9]+\\.([0-9]+).*" "\\1" MUPARSER_MINOR_VERSION ${MUPARSER_VERSION_STRING} )
      string ( REGEX REPLACE "[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" MUPARSER_PATCH_VERSION ${MUPARSER_VERSION_STRING} )
      if ( MUPARSER_PATCH_VERSION STREQUAL MUPARSER_VERSION_STRING )
        set ( MUPARSER_PATCH_VERSION )
      endif ()
    endif ()
  endif ()
endif ()

# check version
set ( _MUPARSER_VERSION_MATCH TRUE )
if ( MuParser_FIND_VERSION AND MUPARSER_VERSION_STRING )
  if ( MuParser_FIND_VERSION_EXACT )
    if ( ${MuParser_FIND_VERSION} VERSION_EQUAL ${MUPARSER_VERSION_STRING} )
    else()
      set ( _MUPARSER_VERSION_MATCH FALSE)
    endif ()
  else ()
    if ( ${MuParser_FIND_VERSION} VERSION_GREATER ${MUPARSER_VERSION_STRING} )
      set ( _MUPARSER_VERSION_MATCH FALSE )
    endif ()
  endif ()
endif ()

find_library ( MUPARSER_LIBRARY
               NAMES muparser
               PATH_SUFFIXES muparser )
               
# set MUPARSER_LIBRARIES
set ( MUPARSER_LIBRARIES ${MUPARSER_LIBRARY} )

# root dir
# try to guess root dir from include dir
if ( MUPARSER_INCLUDE_DIR )
  string ( REGEX REPLACE "(.*)/include.*" "\\1" MUPARSER_ROOT_DIR ${MUPARSER_INCLUDE_DIR} )

# try to guess root dir from library dir
elseif ( MUPARSER_LIBRARY )
  string ( REGEX REPLACE "(.*)/lib[/|32|64].*" "\\1" MUPARSER_ROOT_DIR ${MUPARSER_LIBRARY} )
endif ()

# handle REQUIRED and QUIET options
include ( FindPackageHandleStandardArgs )

if ( CMAKE_VERSION VERSION_LESS 2.8.3 )
  find_package_handle_standard_args ( muParser DEFAULT_MSG MUPARSER_LIBRARY
    _MUPARSER_VERSION_MATCH
    MUPARSER_LIBRARIES
    MUPARSER_INCLUDE_DIR
    MUPARSER_INCLUDE_DIRS
    MUPARSER_VERSION_STRING
    MUPARSER_MAJOR_VERSION
    MUPARSER_MINOR_VERSION
  )
else ()
  find_package_handle_standard_args ( muParser 
    REQUIRED_VARS MUPARSER_LIBRARY _MUPARSER_VERSION_MATCH MUPARSER_LIBRARIES MUPARSER_INCLUDE_DIR MUPARSER_INCLUDE_DIRS
                  MUPARSER_MAJOR_VERSION MUPARSER_MINOR_VERSION
    VERSION_VAR MUPARSER_VERSION_STRING
  )
endif ()

mark_as_advanced (
  MUPARSER_LIBRARY
  MUPARSER_LIBRARIES
  MUPARSER_INCLUDE_DIR
  MUPARSER_INCLUDE_DIRS
  MUPARSER_ROOT_DIR 
  MUPARSER_VERSION_STRING
  MUPARSER_MAJOR_VERSION
  MUPARSER_MINOR_VERSION
  MUPARSER_PATCH_VERSION
)
