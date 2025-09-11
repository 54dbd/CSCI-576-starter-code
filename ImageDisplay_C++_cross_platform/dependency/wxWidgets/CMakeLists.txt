#############################################################################
# Name:        CMakeLists.txt
# Purpose:     Primary CMake for wxWidgets
# Author:      Tobias Taschner
# Created:     2016-09-20
# Copyright:   (c) 2016 wxWidgets development team
# Licence:     wxWindows licence
#############################################################################

cmake_minimum_required(VERSION 4.0)

if(NOT CMAKE_CONFIGURATION_TYPES)
    get_property(HAVE_MULTI_CONFIG_GENERATOR GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    # Set default configuration types for multi-config generators
    if(HAVE_MULTI_CONFIG_GENERATOR)
        set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
    endif()
endif()

file(TO_CMAKE_PATH "${CMAKE_INSTALL_PREFIX}" CMAKE_INSTALL_PREFIX)

# https://blog.kitware.com/cmake-and-the-default-build-type/
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(default_build_type "Debug")
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# This block, particularly the versions used, should be kept in sync with
# samples/minimal/CMakeLists.txt.
if(APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
    # If no deployment target has been set default to the minimum supported
    # OS version (this has to be set before the first project() call)
    if(CMAKE_SYSTEM_NAME STREQUAL "iOS")
        set(CMAKE_OSX_DEPLOYMENT_TARGET 12.0 CACHE STRING "iOS Deployment Target")
    else()
        set(CMAKE_OSX_DEPLOYMENT_TARGET 10.10 CACHE STRING "macOS Deployment Target")
    endif()
endif()

include(build/cmake/policies.cmake NO_POLICY_SCOPE)

# Initialize variables for quick access to wx root dir in sub dirs
set(wxSOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(wxBINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(wxOUTPUT_DIR ${wxBINARY_DIR}/lib)

# parse the version number from wx/version.h and include in wxMAJOR_VERSION and wxMINOR_VERSION
file(READ "${wxSOURCE_DIR}/include/wx/version.h" WX_VERSION_H_CONTENTS)
string(REGEX MATCH "wxMAJOR_VERSION[ \t]+([0-9]+)"
    wxMAJOR_VERSION ${WX_VERSION_H_CONTENTS})
string (REGEX MATCH "([0-9]+)"
    wxMAJOR_VERSION ${wxMAJOR_VERSION})
string(REGEX MATCH "wxMINOR_VERSION[ \t]+([0-9]+)"
    wxMINOR_VERSION ${WX_VERSION_H_CONTENTS})
string (REGEX MATCH "([0-9]+)"
    wxMINOR_VERSION ${wxMINOR_VERSION})
string(REGEX MATCH "wxRELEASE_NUMBER[ \t]+([0-9]+)"
    wxRELEASE_NUMBER ${WX_VERSION_H_CONTENTS})
string (REGEX MATCH "([0-9]+)"
    wxRELEASE_NUMBER ${wxRELEASE_NUMBER})

# Determine if current version is a "Development" release
math(EXPR rel_dev "${wxMINOR_VERSION} % 2")
if(rel_dev)
    set(wxVERSION_IS_DEV TRUE)
else()
    set(wxVERSION_IS_DEV FALSE)
endif()

# parse the .so version from build/bakefiles/version.bkl
file(READ "${wxSOURCE_DIR}/build/bakefiles/version.bkl" WX_VERSION_BKL_CONTENTS)

# --- robust regex: accept any non-digit separator between name and number ---
string(REGEX MATCH "WX_CURRENT[^0-9]*([0-9]+)" _match_current "${WX_VERSION_BKL_CONTENTS}")
if(_match_current)
  set(WX_CURRENT "${CMAKE_MATCH_1}")
else()
  message(WARNING "Failed to parse WX_CURRENT from ${wxSOURCE_DIR}/build/bakefiles/version.bkl. Defaulting WX_CURRENT=0")
  set(WX_CURRENT "0")
endif()

string(REGEX MATCH "WX_REVISION[^0-9]*([0-9]+)" _match_revision "${WX_VERSION_BKL_CONTENTS}")
if(_match_revision)
  set(WX_REVISION "${CMAKE_MATCH_1}")
else()
  message(WARNING "Failed to parse WX_REVISION from ${wxSOURCE_DIR}/build/bakefiles/version.bkl. Defaulting WX_REVISION=0")
  set(WX_REVISION "0")
endif()

string(REGEX MATCH "WX_AGE[^0-9]*([0-9]+)" _match_age "${WX_VERSION_BKL_CONTENTS}")
if(_match_age)
  set(WX_AGE "${CMAKE_MATCH_1}")
else()
  message(WARNING "Failed to parse WX_AGE from ${wxSOURCE_DIR}/build/bakefiles/version.bkl. Defaulting WX_AGE=0")
  set(WX_AGE "0")
endif()

# debug: print parsed values to confirm
message(STATUS "Parsed from version.bkl: WX_CURRENT='${WX_CURRENT}' WX_AGE='${WX_AGE}' WX_REVISION='${WX_REVISION}'")

# Extra safety: strip any non-digit just in case (should be unnecessary after above, but harmless)
string(REGEX REPLACE "[^0-9]" "" WX_CURRENT "${WX_CURRENT}")
string(REGEX REPLACE "[^0-9]" "" WX_AGE "${WX_AGE}")
string(REGEX REPLACE "[^0-9]" "" WX_REVISION "${WX_REVISION}")

# final check before math
if(WX_CURRENT STREQUAL "" OR WX_AGE STREQUAL "")
  message(FATAL_ERROR "WX_CURRENT or WX_AGE is empty after parsing. Can't compute wxSOVERSION_MAJOR. Dumping file contents for debugging:\n${WX_VERSION_BKL_CONTENTS}")
endif()

math(EXPR wxSOVERSION_MAJOR "${WX_CURRENT} - ${WX_AGE}")
set(wxSOVERSION "${wxSOVERSION_MAJOR}.${WX_AGE}.${WX_REVISION}")

message(STATUS "Computed wxSOVERSION = ${wxSOVERSION}")

# preserve your existing version/copyright lines
set(wxVERSION ${wxMAJOR_VERSION}.${wxMINOR_VERSION}.${wxRELEASE_NUMBER})
set(wxCOPYRIGHT "2002-2023 wxWidgets")

set(wxLANGUAGES C CXX)
if(APPLE AND NOT CMAKE_VERSION VERSION_LESS "3.16")
    list(APPEND wxLANGUAGES OBJCXX)
endif()

project(wxWidgets VERSION ${wxVERSION} LANGUAGES ${wxLANGUAGES})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(build/cmake/main.cmake)

# Set the default startup project for Visual Studio
if(wxBUILD_SAMPLES AND wxUSE_GUI)
    set_directory_properties(PROPERTIES VS_STARTUP_PROJECT minimal)
endif()
