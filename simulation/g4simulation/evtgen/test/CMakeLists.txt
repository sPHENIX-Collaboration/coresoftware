
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


# the test executables require ROOT for histogramming etc.
find_package(ROOT QUIET COMPONENTS Hist)
if (ROOT_FOUND)
    add_library(root_interface INTERFACE IMPORTED)
    target_include_directories(root_interface INTERFACE ${ROOT_INCLUDE_DIRS})
    target_link_libraries(root_interface INTERFACE ${ROOT_LIBRARIES})

    # build each of the executables
    foreach( test_exe evt_dalitz evtgenlhc_test1 example1 exampleWriteHepMC )
        add_executable(${test_exe} ${test_exe}.cc)
        target_link_libraries(${test_exe} PRIVATE EvtGen)
        if (EVTGEN_HEPMC3)
            target_compile_definitions(${test_exe} PRIVATE EVTGEN_HEPMC3)
        endif()
        if( EVTGEN_PYTHIA OR EVTGEN_PHOTOS OR EVTGEN_TAUOLA )
            target_compile_definitions(${test_exe} PRIVATE EVTGEN_EXTERNAL)
            target_link_libraries(${test_exe} PRIVATE EvtGenExternal)
        endif()
        target_link_libraries(${test_exe} PRIVATE root_interface)
    endforeach()

    # install the executables
    install(TARGETS evt_dalitz evtgenlhc_test1 example1 exampleWriteHepMC
        DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/test
        )

    # install the decay files, macros, scripts, etc.
    install(DIRECTORY exampleFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/test)
    install(FILES do_tests
        DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/test
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
else()
    message(WARNING "Could not find ROOT, cannot build test executables")
endif()
