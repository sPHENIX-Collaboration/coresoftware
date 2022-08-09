# building tests
#------------------
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
#


include_directories(${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/search/include ${ROOT_INCLUDE_DIR})

#---Add executables------------------------------------------------------------
add_executable(basic_tree.exe ${CMAKE_CURRENT_SOURCE_DIR}/basic_tree.cc)
target_link_libraries(basic_tree.exe ${HEPMC3_LIB}  ${ROOT_LIBRARIES})
set_target_properties(basic_tree.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/$<0:>)
if (USE_INSTALLED_HEPMC3)
install(TARGETS basic_tree.exe DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if (CMAKE_GENERATOR MATCHES "Visual Studio")
 message(STATUS "HepMC3: Visual Studio does not support targets with mixed code, please use e.g. NMake Files.")
else()
if(CMAKE_Fortran_COMPILER)
 add_executable(hepevt_wrapper_example.exe ${CMAKE_CURRENT_SOURCE_DIR}/hepevt_wrapper_example_main.cc ${CMAKE_CURRENT_SOURCE_DIR}/hepevt_wrapper_example_fortran.f)
 target_link_libraries(hepevt_wrapper_example.exe ${HEPMC3_LIB})
 set_target_properties(hepevt_wrapper_example.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/$<0:>)
if (USE_INSTALLED_HEPMC3)
install(TARGETS hepevt_wrapper_example.exe DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
else()
 message(STATUS "HepMC3: No Fortran compiler - dissabling HepEvt example.")
endif()
endif()

add_executable(HepMC2_reader_example.exe ${CMAKE_CURRENT_SOURCE_DIR}/HepMC2_reader_example.cc)
target_link_libraries(HepMC2_reader_example.exe ${HEPMC3_LIB} ${ROOT_LIBRARIES})
set_target_properties(HepMC2_reader_example.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/$<0:>)

if (USE_INSTALLED_HEPMC3)
install(TARGETS HepMC2_reader_example.exe DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

add_executable(HepMC3_fileIO_example.exe ${CMAKE_CURRENT_SOURCE_DIR}/HepMC3_fileIO_example.cc)
target_link_libraries(HepMC3_fileIO_example.exe ${HEPMC3_LIB} ${ROOT_LIBRARIES})
set_target_properties(HepMC3_fileIO_example.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/$<0:>)

if (USE_INSTALLED_HEPMC3)
install(TARGETS HepMC3_fileIO_example.exe DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
