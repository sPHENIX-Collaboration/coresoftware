find_package(Doxygen)
if (DOXYGEN_FOUND)
  find_program(DOXYGEN NAMES doxygen doxygen.exe HINTS ${DOXYGEN_EXECUTABLE})
else()
  find_program(DOXYGEN NAMES doxygen doxygen.exe )
endif()
if (DOXYGEN)
  configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in"
  "${CMAKE_CURRENT_BINARY_DIR}/doxygen/Doxyfile"
  IMMEDIATE @ONLY)
  message(STATUS "HepMC3 doc: doxygen found in ${DOXYGEN}. Generation of documentation is enabled.")
  add_custom_target( doxydoc ALL
  COMMAND ${DOXYGEN}  ${CMAKE_CURRENT_BINARY_DIR}/doxygen/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen
  )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen/html DESTINATION ${CMAKE_INSTALL_DOCDIR} COMPONENT doc)
else()
  message(STATUS "HepMC3 doc: doxygen not found, generation of documentation is disabled.")
endif()
