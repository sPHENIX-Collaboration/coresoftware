
########################################################################
# Copyright 1998-2020 CERN for the benefit of the EvtGen authors       #
#                                                                      #
# This file is part of EvtGen.                                         #
#                                                                      #
# EvtGen is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# EvtGen is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     #
########################################################################


# Set the build type (if not already specified)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "EvtGen: Setting build type to 'Release' as none was specified")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Release, MinSizeRel, Debug, RelWithDebInfo" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
elseif(CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "EvtGen: Build type '${CMAKE_BUILD_TYPE}'")
endif()

# Set the warning/optimise/debug flags for each build type
if( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" )
    message(STATUS "EvtGen: Customising compiler flags for each build type")

    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsigned-char -Wall -Wextra -Wshadow -Woverloaded-virtual")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsigned-char -Wall -Wextra")

    if( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" )
        set(CMAKE_CXX_FLAGS_DEBUG          "-Og -g3")
        set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
        set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g3")
    elseif( ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" )
        set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
        set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
        set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    endif()
else()
    message(STATUS "EvtGen: No customisation of compiler flags implemented for compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

# Make sure our project's include directories always come first
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)

# Prioritise UNIX-style packages over macOS frameworks
set(CMAKE_FIND_FRAMEWORK LAST)

# Control verbosity of the build
set(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Control verbosity of generated Makefiles")

# C++ standard settings
set(CMAKE_CXX_EXTENSIONS OFF)
if(DEFINED ENV{CMAKE_CXX_STANDARD})
    set(CMAKE_CXX_STANDARD $ENV{CMAKE_CXX_STANDARD} CACHE STRING "C++ standard")
else()
    set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ standard")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "EvtGen: Using C++${CMAKE_CXX_STANDARD} standard")

# Special linker flags for MacOSX
if (APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -single_module -undefined dynamic_lookup")
endif()

# RPATH handling
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

