
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

# - Try to find Tauola++
# Defines:
#
#  Tauola++_FOUND
#  Tauola++_INCLUDE_DIR
#  Tauola++_INCLUDE_DIRS (not cached)
#  Tauola++_<component>_LIBRARY
#  Tauola++_<component>_FOUND
#  Tauola++_LIBRARIES (not cached)
#  Tauola++_LIBRARY_DIRS (not cached)

# Check all if none is explicitly requested
if(NOT Tauola++_FIND_COMPONENTS)
  set(Tauola++_FIND_COMPONENTS Fortran CxxInterface HepMC HepMC3)
endif()
set(Tauola++_FOUND TRUE)
foreach(component ${Tauola++_FIND_COMPONENTS})
  find_library(Tauola++_${component}_LIBRARY NAMES Tauola${component}
               HINTS ${Tauola++_ROOT_DIR}/lib $ENV{TAUOLAPP_ROOT_DIR}/lib ${TAUOLAPP_ROOT_DIR}/lib
                     ${Tauola++_ROOT_DIR}/lib64 $ENV{TAUOLAPP_ROOT_DIR}/lib64 ${TAUOLAPP_ROOT_DIR}/lib64
                      )
  if (Tauola++_${component}_LIBRARY)
    set(Tauola++_${component}_FOUND TRUE)
    list(APPEND Tauola++_LIBRARIES ${Tauola++_${component}_LIBRARY})
    get_filename_component(libdir ${Tauola++_${component}_LIBRARY} PATH)
    list(APPEND Tauola++_LIBRARY_DIRS ${libdir})
  else()
    set(Tauola++_${component}_FOUND FALSE)
    if (Tauola++_FIND_REQUIRED_${component})
    set(Tauola++_FOUND FALSE)
    endif()
  endif()
  message(STATUS "EvtGen: Tauola++_${component}_FOUND=${Tauola++_${component}_FOUND}, Tauola++_FIND_REQUIRED_${component}=${Tauola++_FIND_REQUIRED_${component}}  in ${Tauola++_${component}_LIBRARY}")
  mark_as_advanced(Tauola++_${component}_LIBRARY)
endforeach()

if(Tauola++_LIBRARY_DIRS)
  list(REMOVE_DUPLICATES Tauola++_LIBRARY_DIRS)
endif()

find_path(Tauola++_INCLUDE_DIR Tauola/Tauola.h
          HINTS ${Tauola++_ROOT_DIR}/include
                $ENV{TAUOLAPP_ROOT_DIR}/include ${TAUOLAPP_ROOT_DIR}/include)
set(Tauola++_INCLUDE_DIRS ${Tauola++_INCLUDE_DIR})
mark_as_advanced(Tauola++_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set Tauola++_FOUND to TRUE if
# all listed variables are TRUE
if (Tauola++_FOUND)
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Tauola++ DEFAULT_MSG Tauola++_INCLUDE_DIR Tauola++_LIBRARIES)
endif()
mark_as_advanced(Tauola++_FOUND)
