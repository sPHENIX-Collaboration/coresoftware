
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

# Top level CMakeLists.txt for EvtGen

# Enforce an out-of-source build.
# Should be the first action in the top level CMakeLists.txt
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "EvtGen requires an out-of-source build.")
    message(STATUS "Please remove the following files from ${CMAKE_BINARY_DIR}:")
    message(STATUS "    CMakeCache.txt")
    message(STATUS "    CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there, pointing it to:")
    message(STATUS "    ${CMAKE_SOURCE_DIR}")
    message(FATAL_ERROR "in-source build detected")
endif()

# Also require a minimum version of CMake
cmake_minimum_required(VERSION 3.11.0)

# Project setup
if(${CMAKE_VERSION} VERSION_LESS 3.12.0)
    project(EvtGen
        VERSION 2.1.1
        DESCRIPTION "Monte Carlo generator of particle decays, in particular the weak decays of heavy flavour particles such as B mesons."
        LANGUAGES C CXX
        )
else()
    project(EvtGen
        VERSION 2.1.1
        DESCRIPTION "Monte Carlo generator of particle decays, in particular the weak decays of heavy flavour particles such as B mesons."
        LANGUAGES C CXX
        HOMEPAGE_URL "https://evtgen.hepforge.org"
        )
endif()

# Prepend this project's custom module path(s) to CMAKE_MODULE_PATH
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
# Include needed modules to perform any custom setup
# Install paths
include(GNUInstallDirs)
# Compilation/linking flags and related settings
include(CompilerFlags)
# EvtGen external dependencies
option(EVTGEN_HEPMC3 "Enable/disable linking with HepMC3"   ON)
option(EVTGEN_PYTHIA "Enable/disable linking with Pythia8"  OFF)
option(EVTGEN_PHOTOS "Enable/disable linking with Photos++" OFF)
option(EVTGEN_TAUOLA "Enable/disable linking with Tauola++" OFF)
option(EVTGEN_SUPPRESS_EXTERNAL_WARNINGS "Enable/disable suppression of warnings from external dependencies" OFF)
if (EVTGEN_HEPMC3)
    message(STATUS "EvtGen: Linking with HepMC3")
else()
    message(STATUS "EvtGen: Linking with HepMC2")
endif()
message(STATUS "EvtGen: Optional linking with HepMC3   EVTGEN_HEPMC3 ${EVTGEN_HEPMC3}")
message(STATUS "EvtGen: Optional linking with Pythia8  EVTGEN_PYTHIA ${EVTGEN_PYTHIA}")
message(STATUS "EvtGen: Optional linking with Photos++ EVTGEN_PHOTOS ${EVTGEN_PHOTOS}")
message(STATUS "EvtGen: Optional linking with Tauola++ EVTGEN_TAUOLA ${EVTGEN_TAUOLA}")
include(ExternalDependencies)

# Now build the library
add_subdirectory(src)

# Copy the particle property and decay tables
configure_file(DECAY.DEC ${CMAKE_CURRENT_BINARY_DIR}/DECAY.DEC COPYONLY)
configure_file(DECAY.XML ${CMAKE_CURRENT_BINARY_DIR}/DECAY.XML COPYONLY)
configure_file(evt.pdl ${CMAKE_CURRENT_BINARY_DIR}/evt.pdl COPYONLY)
# Copy the older decay tables
configure_file(DECAY_2009.DEC ${CMAKE_CURRENT_BINARY_DIR}/DECAY_2009.DEC COPYONLY)
configure_file(DECAY_2009.XML ${CMAKE_CURRENT_BINARY_DIR}/DECAY_2009.XML COPYONLY)
configure_file(DECAY_2010.DEC ${CMAKE_CURRENT_BINARY_DIR}/DECAY_2010.DEC COPYONLY)
configure_file(DECAY_2010.XML ${CMAKE_CURRENT_BINARY_DIR}/DECAY_2010.XML COPYONLY)

# Build the executables in the test and validation directories
option(EVTGEN_BUILD_DOC         "Enable/disable building of userguide documentation in 'doc' directory"   OFF)
option(EVTGEN_BUILD_DOXYGEN     "Enable/disable building of Doxygen documentation in 'doxygen' directory" OFF)
option(EVTGEN_BUILD_TESTS       "Enable/disable building of executables in 'test' directory"              OFF)
option(EVTGEN_BUILD_VALIDATIONS "Enable/disable building of executables in 'validation' directory"        OFF)
message(STATUS "EvtGen: Building of userguide documentation in 'doc' directory   EVTGEN_BUILD_DOC         ${EVTGEN_BUILD_DOC}")
message(STATUS "EvtGen: Building of Doxygen documentation in 'doxygen' directory EVTGEN_BUILD_DOXYGEN     ${EVTGEN_BUILD_DOXYGEN}")
message(STATUS "EvtGen: Building of executables in 'test' directory              EVTGEN_BUILD_TESTS       ${EVTGEN_BUILD_TESTS}")
message(STATUS "EvtGen: Building of executables in 'validation' directory        EVTGEN_BUILD_VALIDATIONS ${EVTGEN_BUILD_VALIDATIONS}")
if(EVTGEN_BUILD_DOC)
    add_subdirectory(doc)
endif()
if(EVTGEN_BUILD_DOXYGEN)
    add_subdirectory(doxygen)
endif()
if(EVTGEN_BUILD_TESTS)
    add_subdirectory(test)
endif()
if(EVTGEN_BUILD_VALIDATIONS)
    add_subdirectory(validation)
endif()

# Install the include directories
install(DIRECTORY EvtGen         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY EvtGenBase     DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY EvtGenExternal DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY EvtGenModels   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install the particle properties and decay tables
install(FILES DECAY.DEC DECAY.XML evt.pdl DECAY_2010.XML DECAY_2010.DEC DECAY_2009.XML DECAY_2009.DEC DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME})

# Install the COPYING and AUTHORS files
install(FILES AUTHORS COPYING DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME})

# Install the README and History files
install(FILES README.md Pythia8_README.md Tauola_README.md History.md DESTINATION ${CMAKE_INSTALL_DOCDIR})

# Generate CMake config files, which can be used by other projects
include(CMakePackageConfigHelpers)

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
set(DATA_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME})

configure_package_config_file(cmake/Templates/EvtGenConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/EvtGenConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/cmake
    PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR DATA_INSTALL_DIR
    )

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/EvtGenConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion
    )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/EvtGenConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/EvtGenConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/cmake
    )

install(
    EXPORT "EvtGenTargets"
    NAMESPACE "EvtGen::"
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/cmake
    )
