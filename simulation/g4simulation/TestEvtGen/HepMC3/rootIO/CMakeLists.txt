#--------
#-root IO

include_directories(${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Generate dictionaries
#
if (MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR})
else()
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR})
endif()
ROOT_GENERATE_DICTIONARY(G__HepMC3rootIO include/rootIO_Classes.hh LINKDEF include/LinkDef.hh NOINSTALL OPTIONS -inlineInputHeader)

add_library(HepMC3rootIO SHARED ${PROJECT_SOURCE_DIR}/rootIO/src/WriterRoot.cc ${PROJECT_SOURCE_DIR}/rootIO/src/ReaderRoot.cc ${PROJECT_SOURCE_DIR}/rootIO/src/WriterRootTree.cc ${PROJECT_SOURCE_DIR}/rootIO/src/ReaderRootTree.cc ${PROJECT_SOURCE_DIR}/rootIO/src/Streamers.cc G__HepMC3rootIO.cxx)
set_property(TARGET HepMC3rootIO PROPERTY POSITION_INDEPENDENT_CODE 1)

target_link_libraries(HepMC3rootIO ${ROOT_Tree_LIBRARY} ${ROOT_RIO_LIBRARY} ${ROOT_Core_LIBRARY} HepMC3)
set_target_properties(HepMC3rootIO PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>)
if (MSVC)
  set_target_properties(HepMC3rootIO PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>)
endif()
set_target_properties(HepMC3rootIO PROPERTIES SOVERSION 3)
# installs
install(TARGETS HepMC3rootIO DESTINATION ${HEPMC3_ROOTIO_INSTALL_LIBDIR} COMPONENT rootIOlibs)
install(DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/ DESTINATION  ${HEPMC3_ROOTIO_INSTALL_LIBDIR}  COMPONENT rootIOlibs  FILES_MATCHING  PATTERN "*.pcm"  PATTERN "*.rootmap"  PATTERN "CMakeFiles" EXCLUDE   PATTERN "src" EXCLUDE  PATTERN "include" EXCLUDE)
install(DIRECTORY include/HepMC3 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT rootIOdevel)
