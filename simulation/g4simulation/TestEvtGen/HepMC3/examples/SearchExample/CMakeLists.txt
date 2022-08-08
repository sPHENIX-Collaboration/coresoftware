# building tests
#------------------
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
#

#---Add executables------------------------------------------------------------
add_executable(search_example.exe ${CMAKE_CURRENT_SOURCE_DIR}/src/search_example.cc)
target_link_libraries(search_example.exe ${HEPMC3_LIB} ${HEPMC3_SEARCH_LIB} )
target_include_directories(search_example.exe PUBLIC  ${CMAKE_SOURCE_DIR}/search/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(search_example.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/$<0:>)
if (USE_INSTALLED_HEPMC3)
  install(TARGETS search_example.exe DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
