
include_directories(${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}  ${PROJECT_SOURCE_DIR}/rootIO/include ${ROOT_INCLUDE_DIR})

add_executable(rootIOTree_example_read.exe ${CMAKE_CURRENT_SOURCE_DIR}/rootIOTree_example_read.cc)
target_link_libraries(rootIOTree_example_read.exe ${ROOT_LIBRARIES}  ${HEPMC3_LIB} ${HEPMC3_ROOTIO_LIB} )
set_target_properties(rootIOTree_example_read.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR})

add_executable(rootIOTree_example_write.exe ${CMAKE_CURRENT_SOURCE_DIR}/rootIOTree_example_write.cc)
target_link_libraries(rootIOTree_example_write.exe ${ROOT_LIBRARIES}  ${HEPMC3_LIB} ${HEPMC3_ROOTIO_LIB} )
set_target_properties(rootIOTree_example_write.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR})

# create environment scripts

set(bindir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")
set(libdir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:${ROOT_LIBRARY_DIR}")

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/rootIOTree_example_env.sh.in
                ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/rootIOTree_example_env.sh
                @ONLY
              )

# installs
if (USE_INSTALLED_HEPMC3)
  install(TARGETS rootIOTree_example_write.exe rootIOTree_example_read.exe DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/ DESTINATION ${CMAKE_INSTALL_BINDIR}
        FILES_MATCHING PATTERN "*env.sh*")
endif()
