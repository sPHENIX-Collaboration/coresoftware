
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


# Use glob to find the sources for the main and external libraries
file(GLOB EVTGEN_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/EvtGenBase/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/EvtGenModels/*.cpp
    )

file(GLOB EVTGEN_EXTERNAL_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/EvtGenExternal/*.cpp
    )

# Add the main EvtGen library...
# First build all the object files (with PIC enabled) so they can be used to
# build both the shared and static libs
add_library(objlib OBJECT ${EVTGEN_SOURCES})
set_target_properties(objlib PROPERTIES POSITION_INDEPENDENT_CODE 1)
target_compile_definitions(objlib PRIVATE EVTGEN_CPP11)
target_include_directories(objlib PRIVATE ${CMAKE_SOURCE_DIR})
if(EVTGEN_HEPMC3)
    target_compile_definitions(objlib PRIVATE EVTGEN_HEPMC3)
    if (EVTGEN_SUPPRESS_EXTERNAL_WARNINGS)
        target_include_directories(objlib SYSTEM PRIVATE ${HEPMC3_INCLUDE_DIR})
    else()
        target_include_directories(objlib PRIVATE ${HEPMC3_INCLUDE_DIR})
    endif()
else()
    if (EVTGEN_SUPPRESS_EXTERNAL_WARNINGS)
        target_include_directories(objlib SYSTEM PRIVATE ${HEPMC2_INCLUDE_DIR})
    else()
        target_include_directories(objlib PRIVATE ${HEPMC2_INCLUDE_DIR})
    endif()
endif()

# Now make the shared library from the object files
add_library(EvtGen SHARED $<TARGET_OBJECTS:objlib>)
set_target_properties(EvtGen PROPERTIES OUTPUT_NAME EvtGen)
set_target_properties(EvtGen PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR} )
set_target_properties(EvtGen PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
target_compile_definitions(EvtGen PUBLIC EVTGEN_CPP11)
target_include_directories(EvtGen PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if(EVTGEN_HEPMC3)
    target_compile_definitions(EvtGen PUBLIC EVTGEN_HEPMC3)
    target_include_directories(EvtGen PUBLIC ${HEPMC3_INCLUDE_DIR})
    target_link_libraries(EvtGen PUBLIC ${HEPMC3_LIB} ${HEPMC3_SEARCH_LIB})
else()
    target_include_directories(EvtGen PUBLIC ${HEPMC2_INCLUDE_DIR})
    target_link_libraries(EvtGen PUBLIC ${HEPMC2_LIBRARIES})
endif()

# Now make the static library from the object files
add_library(EvtGenStatic STATIC $<TARGET_OBJECTS:objlib>)
set_target_properties(EvtGenStatic PROPERTIES OUTPUT_NAME EvtGen)
set_target_properties(EvtGenStatic PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
target_compile_definitions(EvtGenStatic PUBLIC EVTGEN_CPP11)
target_include_directories(EvtGenStatic PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if(EVTGEN_HEPMC3)
    target_compile_definitions(EvtGenStatic PUBLIC EVTGEN_HEPMC3)
    target_include_directories(EvtGenStatic PUBLIC ${HEPMC3_INCLUDE_DIR})
    target_link_libraries(EvtGenStatic PUBLIC ${HEPMC3_LIB} ${HEPMC3_SEARCH_LIB})
else()
    target_include_directories(EvtGenStatic PUBLIC ${HEPMC2_INCLUDE_DIR})
    target_link_libraries(EvtGenStatic PUBLIC ${HEPMC2_LIBRARIES})
endif()



# Add the EvtGenExternal library (if required)...
if( EVTGEN_PYTHIA OR EVTGEN_PHOTOS OR EVTGEN_TAUOLA )
    # Again, first build all the object files (with PIC enabled) so they can be
    # used to build both the shared and static libs
    add_library(objlib_ext OBJECT ${EVTGEN_EXTERNAL_SOURCES})
    set_target_properties(objlib_ext PROPERTIES POSITION_INDEPENDENT_CODE 1)
    target_compile_definitions(objlib_ext PRIVATE EVTGEN_CPP11)
    target_include_directories(objlib_ext PRIVATE ${CMAKE_SOURCE_DIR})
    if(EVTGEN_HEPMC3)
        target_compile_definitions(objlib_ext PRIVATE EVTGEN_HEPMC3)
        if (EVTGEN_SUPPRESS_EXTERNAL_WARNINGS)
            target_include_directories(objlib_ext SYSTEM PRIVATE ${HEPMC3_INCLUDE_DIR})
        else()
            target_include_directories(objlib_ext PRIVATE ${HEPMC3_INCLUDE_DIR})
        endif()
    else()
        if (EVTGEN_SUPPRESS_EXTERNAL_WARNINGS)
            target_include_directories(objlib_ext SYSTEM PRIVATE ${HEPMC2_INCLUDE_DIR})
        else()
            target_include_directories(objlib_ext PRIVATE ${HEPMC2_INCLUDE_DIR})
        endif()
    endif()
    if(EVTGEN_PYTHIA)
        target_compile_definitions(objlib_ext PRIVATE EVTGEN_PYTHIA)
        if (EVTGEN_SUPPRESS_EXTERNAL_WARNINGS)
            target_include_directories(objlib_ext SYSTEM PRIVATE ${PYTHIA8_INCLUDE_DIRS})
        else()
            target_include_directories(objlib_ext PRIVATE ${PYTHIA8_INCLUDE_DIRS})
        endif()
    endif()
    if(EVTGEN_PHOTOS)
        target_compile_definitions(objlib_ext PRIVATE EVTGEN_PHOTOS)
        if (EVTGEN_SUPPRESS_EXTERNAL_WARNINGS)
            target_include_directories(objlib_ext SYSTEM PRIVATE ${Photos++_INCLUDE_DIRS})
        else()
            target_include_directories(objlib_ext PRIVATE ${Photos++_INCLUDE_DIRS})
        endif()
    endif()
    if(EVTGEN_TAUOLA)
        target_compile_definitions(objlib_ext PRIVATE EVTGEN_TAUOLA)
        if (EVTGEN_SUPPRESS_EXTERNAL_WARNINGS)
            target_include_directories(objlib_ext SYSTEM PRIVATE ${Tauola++_INCLUDE_DIRS})
        else()
            target_include_directories(objlib_ext PRIVATE ${Tauola++_INCLUDE_DIRS})
        endif()
    endif()

    # Now make the shared library from the object files
    add_library(EvtGenExternal SHARED $<TARGET_OBJECTS:objlib_ext>)
    set_target_properties(EvtGenExternal PROPERTIES OUTPUT_NAME EvtGenExternal)
    set_target_properties(EvtGenExternal PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR} )
    set_target_properties(EvtGenExternal PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
    target_compile_definitions(EvtGenExternal PUBLIC EVTGEN_CPP11)
    target_include_directories(EvtGenExternal PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    if(EVTGEN_PYTHIA)
        target_compile_definitions(EvtGenExternal PUBLIC EVTGEN_PYTHIA)
        target_include_directories(EvtGenExternal PUBLIC ${PYTHIA8_INCLUDE_DIRS})
        target_link_libraries(EvtGenExternal PUBLIC ${PYTHIA8_LIBRARIES})
    endif()
    if(EVTGEN_HEPMC3)
        target_compile_definitions(EvtGenExternal PUBLIC EVTGEN_HEPMC3)
        target_include_directories(EvtGenExternal PUBLIC ${HEPMC3_INCLUDE_DIR})
        target_link_libraries(EvtGenExternal PUBLIC ${HEPMC3_LIB} ${HEPMC3_SEARCH_LIB})
        if(EVTGEN_PHOTOS)
            target_compile_definitions(EvtGenExternal PUBLIC EVTGEN_PHOTOS)
            target_include_directories(EvtGenExternal PUBLIC ${Photos++_INCLUDE_DIRS})
            # From version 3.64 Photos has HepMC3 support
            target_link_libraries(EvtGenExternal PUBLIC ${Photos++_pp_LIBRARY} ${Photos++_ppHepMC3_LIBRARY})
        endif()
        if(EVTGEN_TAUOLA)
            target_compile_definitions(EvtGenExternal PUBLIC EVTGEN_TAUOLA)
            target_include_directories(EvtGenExternal PUBLIC ${Tauola++_INCLUDE_DIRS})
            # From version 1.1.8 Tauola has HepMC3 support
            target_link_libraries(EvtGenExternal PUBLIC ${Tauola++_CxxInterface_LIBRARY} ${Tauola++_Fortran_LIBRARY} ${Tauola++_HepMC3_LIBRARY})
        endif()
    else()
        target_include_directories(EvtGenExternal PUBLIC ${HEPMC2_INCLUDE_DIR})
        target_link_libraries(EvtGenExternal PUBLIC ${HEPMC2_LIBRARIES})
        if(EVTGEN_PHOTOS)
            target_compile_definitions(EvtGenExternal PUBLIC EVTGEN_PHOTOS)
            target_include_directories(EvtGenExternal PUBLIC ${Photos++_INCLUDE_DIRS})
            # Photos has different library structures for versions before and after 3.58
            # so we need to check which one we have
            if(Photos++_pp_FOUND AND Photos++_ppHepMC_FOUND)
                message(STATUS "EvtGen: PHOTOS has pp and ppHepMC components")
                target_link_libraries(EvtGenExternal PUBLIC ${Photos++_pp_LIBRARY} ${Photos++_ppHepMC_LIBRARY})
            elseif(Photos++_CxxInterface_FOUND AND Photos++_Fortran_FOUND)
                message(STATUS "EvtGen: PHOTOS has CxxInterface and Fortran components")
                target_link_libraries(EvtGenExternal PUBLIC ${Photos++_CxxInterface_LIBRARY} ${Photos++_Fortran_LIBRARY})
            else()
                message(ERROR "EvtGen: PHOTOS lacks pp+ppHepMC or CxxInterface+Fortran components")
            endif()
        endif()
        if(EVTGEN_TAUOLA)
            target_compile_definitions(EvtGenExternal PUBLIC EVTGEN_TAUOLA)
            target_include_directories(EvtGenExternal PUBLIC ${Tauola++_INCLUDE_DIRS})
            target_link_libraries(EvtGenExternal PUBLIC ${Tauola++_CxxInterface_LIBRARY} ${Tauola++_Fortran_LIBRARY})
            # Older versions of Tauola don't have the HepMC component, the HepMC2 interface is in CxxInterface
            if(Tauola++_HepMC_FOUND)
                target_link_libraries(EvtGenExternal PUBLIC ${Tauola++_HepMC_LIBRARY})
            endif()
        endif()
    endif()

    # Now make the static library from the object files
    add_library(EvtGenExternalStatic STATIC $<TARGET_OBJECTS:objlib_ext>)
    set_target_properties(EvtGenExternalStatic PROPERTIES OUTPUT_NAME EvtGenExternal)
    set_target_properties(EvtGenExternalStatic PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
    target_compile_definitions(EvtGenExternalStatic PUBLIC EVTGEN_CPP11)
    target_include_directories(EvtGenExternalStatic PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    if(EVTGEN_PYTHIA)
        target_compile_definitions(EvtGenExternalStatic PUBLIC EVTGEN_PYTHIA)
        target_include_directories(EvtGenExternalStatic PUBLIC ${PYTHIA8_INCLUDE_DIRS})
        target_link_libraries(EvtGenExternalStatic PUBLIC ${PYTHIA8_LIBRARIES})
    endif()
    if(EVTGEN_HEPMC3)
        target_compile_definitions(EvtGenExternalStatic PUBLIC EVTGEN_HEPMC3)
        target_include_directories(EvtGenExternalStatic PUBLIC ${HEPMC3_INCLUDE_DIR})
        target_link_libraries(EvtGenExternalStatic PUBLIC ${HEPMC3_LIB} ${HEPMC3_SEARCH_LIB})
        if(EVTGEN_PHOTOS)
            target_compile_definitions(EvtGenExternalStatic PUBLIC EVTGEN_PHOTOS)
            target_include_directories(EvtGenExternalStatic PUBLIC ${Photos++_INCLUDE_DIRS})
            # From version 3.64 Photos has HepMC3 support
            target_link_libraries(EvtGenExternalStatic PUBLIC ${Photos++_pp_LIBRARY} ${Photos++_ppHepMC3_LIBRARY})
        endif()
        if(EVTGEN_TAUOLA)
            target_compile_definitions(EvtGenExternalStatic PUBLIC EVTGEN_TAUOLA)
            target_include_directories(EvtGenExternalStatic PUBLIC ${Tauola++_INCLUDE_DIRS})
            # From version 1.1.8 Tauola has HepMC3 support
            target_link_libraries(EvtGenExternalStatic PUBLIC ${Tauola++_CxxInterface_LIBRARY} ${Tauola++_Fortran_LIBRARY} ${Tauola++_HepMC3_LIBRARY})
        endif()
    else()
        target_include_directories(EvtGenExternalStatic PUBLIC ${HEPMC2_INCLUDE_DIR})
        target_link_libraries(EvtGenExternalStatic PUBLIC ${HEPMC2_LIBRARIES})
        if(EVTGEN_PHOTOS)
            target_compile_definitions(EvtGenExternalStatic PUBLIC EVTGEN_PHOTOS)
            target_include_directories(EvtGenExternalStatic PUBLIC ${Photos++_INCLUDE_DIRS})
            # Photos has different library structures for versions before and after 3.58
            # so we need to check which one we have
            if(Photos++_pp_FOUND AND Photos++_ppHepMC_FOUND)
                #message(STATUS "EvtGen: PHOTOS has pp and ppHepMC components")
                target_link_libraries(EvtGenExternalStatic PUBLIC ${Photos++_pp_LIBRARY} ${Photos++_ppHepMC_LIBRARY})
            elseif(Photos++_CxxInterface_FOUND AND Photos++_Fortran_FOUND)
                #message(STATUS "EvtGen: PHOTOS has CxxInterface and Fortran components")
                target_link_libraries(EvtGenExternalStatic PUBLIC ${Photos++_CxxInterface_LIBRARY} ${Photos++_Fortran_LIBRARY})
            else()
                message(ERROR "EvtGen: PHOTOS lacks pp+ppHepMC or CxxInterface+Fortran components")
            endif()
        endif()
        if(EVTGEN_TAUOLA)
            target_compile_definitions(EvtGenExternalStatic PUBLIC EVTGEN_TAUOLA)
            target_include_directories(EvtGenExternalStatic PUBLIC ${Tauola++_INCLUDE_DIRS})
            target_link_libraries(EvtGenExternalStatic PUBLIC ${Tauola++_CxxInterface_LIBRARY} ${Tauola++_Fortran_LIBRARY})
            # Older versions of Tauola don't have the HepMC component, the HepMC2 interface is in CxxInterface
            if(Tauola++_HepMC_FOUND)
                target_link_libraries(EvtGenExternalStatic PUBLIC ${Tauola++_HepMC_LIBRARY})
            endif()
        endif()
    endif()
endif()

# Install the libraries
install(
    TARGETS EvtGen EvtGenStatic
    EXPORT "EvtGenTargets"
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

if( EVTGEN_PYTHIA OR EVTGEN_PHOTOS OR EVTGEN_TAUOLA )
    install(
        TARGETS EvtGenExternal EvtGenExternalStatic
        EXPORT "EvtGenTargets"
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
endif()
