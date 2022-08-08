
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


# the validation executables require ROOT for histogramming etc.
find_package(ROOT QUIET COMPONENTS Hist)
if (ROOT_FOUND)
    add_library(root_interface INTERFACE IMPORTED)
    target_include_directories(root_interface INTERFACE ${ROOT_INCLUDE_DIRS})
    target_link_libraries(root_interface INTERFACE ${ROOT_LIBRARIES})

    # build each of the executables
    foreach( validation_exe compareRootFiles genExampleRootFiles genRootDecayChain testCPVDecays )
        add_executable(${validation_exe} ${validation_exe}.cc)
        target_link_libraries(${validation_exe} PRIVATE EvtGen)
        if (EVTGEN_HEPMC3)
            target_compile_definitions(${validation_exe} PRIVATE EVTGEN_HEPMC3)
        endif()
        if( EVTGEN_PYTHIA OR EVTGEN_PHOTOS OR EVTGEN_TAUOLA )
            target_compile_definitions(${validation_exe} PRIVATE EVTGEN_EXTERNAL)
            target_link_libraries(${validation_exe} PRIVATE EvtGenExternal)
        endif()
        target_link_libraries(${validation_exe} PRIVATE root_interface)
    endforeach()

    # install the executables
    install(TARGETS compareRootFiles genExampleRootFiles genRootDecayChain testCPVDecays
        RUNTIME DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation
        )

    # install the decay files, macros, scripts, etc.
    install(DIRECTORY B0Files DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(DIRECTORY CPVDecayFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(DIRECTORY UpsilonFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(DIRECTORY TauolaFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(DIRECTORY Bs0Files DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(DIRECTORY tauFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(DIRECTORY rootFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(DIRECTORY PHSPFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(DIRECTORY BpFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(DIRECTORY gifFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(DIRECTORY DalitzFiles DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(FILES BKstarGamma.dec BuDst0rhop.dec Kspimumu.dec photosTest.dec DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(FILES compareDalitz.C PhaseSpacePlots.C photosPlots.C plotBKstarGamma.C PlotKspimumu.C DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation)
    install(FILES compareAllDecays.sh
        compareB0PythiaDecays.sh
        compareBpPythiaDecays.sh
        compareBs0PythiaDecays.sh
        compareTauolaDecays.sh
        compareTauPythiaDecays.sh
        compareUpsilonPythiaDecays.sh
        genAllDecayExamples.sh
        genB0PythiaDecays.sh
        genBpPythiaDecays.sh
        genBs0PythiaDecays.sh
        genDalitzDecays.sh
        genDDalitzModes.sh
        genPHSP.sh
        genTauolaDecays.sh
        genTauPythiaDecays.sh
        genUpsilonPythiaDecays.sh
        runKspimumu.sh
        runPhotosTest.sh
        DESTINATION ${CMAKE_INSTALL_DATADIR}/EvtGen/validation
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
else()
    message(WARNING "Could not find ROOT, cannot build validation executables")
endif()
