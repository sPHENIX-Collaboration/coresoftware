
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

# - Locate HepMC library
# in a directory defined via HEPMC2_ROOT_DIR or HEPMC2_DIR environment variable
# Defines:
#
#  HEPMC2_FOUND
#  HEPMC2_INCLUDE_DIR
#  HEPMC2_INCLUDE_DIRS (not cached)
#  HEPMC2_LIBRARIES
#  HEPMC2_FIO_LIBRARIES
set(TEST_HEPMC2_ROOT_DIR  "" ${HEPMC2_ROOT_DIR})
IF(TEST_HEPMC2_ROOT_DIR STREQUAL "")
IF(DEFINED ENV{HEPMC2_ROOT_DIR})
set(HEPMC2_ROOT_DIR  $ENV{HEPMC2_ROOT_DIR})
else()
set(HEPMC2_ROOT_DIR  "/usr")
endif()
endif()

find_path(HEPMC2_INCLUDE_DIR HepMC/GenEvent.h
          HINTS $ENV{HEPMC2_ROOT_DIR}/include ${HEPMC2_ROOT_DIR}/include
          $ENV{HEPMC2_DIR}/include ${HEPMC2_DIR}/include)

find_library(HEPMC2_LIBRARIES NAMES HepMC
             HINTS $ENV{HEPMC2_ROOT_DIR}/lib ${HEPMC2_ROOT_DIR}/lib
             HINTS $ENV{HEPMC2_DIR}/lib ${HEPMC2_DIR}/lib
             HINTS $ENV{HEPMC2_ROOT_DIR}/lib64 ${HEPMC2_ROOT_DIR}/lib64
             HINTS $ENV{HEPMC2_DIR}/lib64 ${HEPMC2_DIR}/lib64
             )

get_filename_component(HEPMC2_LIBRARY_DIR ${HEPMC2_LIBRARIES} PATH)
set(HEPMC2_FIO_LIBRARIES "-L${HEPMC2_LIBRARY_DIR} -lHepMCfio")

set(HEPMC2_INCLUDE_DIRS ${HEPMC2_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set HEPMC2_FOUND to TRUE if
# all listed variables are TRUE

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(HepMC2 FOUND_VAR HEPMC2_FOUND REQUIRED_VARS HEPMC2_INCLUDE_DIR HEPMC2_LIBRARIES)

mark_as_advanced(HEPMC2_FOUND HEPMC2_INCLUDE_DIRS HEPMC2_LIBRARIES)
