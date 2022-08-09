if (USE_INSTALLED_HEPMC3)
#This version is for standalone compilation of examples.
  set(CMAKE_MODULE_PATH    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules    ${CMAKE_MODULE_PATH})
  find_package(Pythia8 8.220)
  if(Pythia8_FOUND)
    message(STATUS "HepMC3 Pythia8Example: Pythia8 package found. pythia8_example enabled")
    add_executable(pythia8_example.exe ${CMAKE_CURRENT_SOURCE_DIR}/pythia8_example.cc )
    if (PYTHIA8_VERSION VERSION_GREATER "8.300")
      target_compile_definitions(pythia8_example.exe PRIVATE "-DHEPMC3_USE_INTERFACE_FROM_PYTHIA8")
      target_include_directories(pythia8_example.exe PRIVATE ${PYTHIA8_INCLUDE_DIR})
    else()
      target_include_directories(pythia8_example.exe PRIVATE ${HEPMC3_INTERFACES_DIR}/pythia8/include/Pythia8 ${PYTHIA8_INCLUDE_DIR})
    endif()
  endif()
else()
#This version is for compilation during installation. Should be consistent with Pythia used in tests
  set(PYTHIA8_ROOT_DIR "" CACHE PATH "Location of Pythia8 installation.")
  find_package(Pythia8 8.220)
  if(Pythia8_FOUND)
    message(STATUS "HepMC3 Pythia8Example: Pythia8 package found. pythia8_example enabled")
    add_executable(pythia8_example.exe ${CMAKE_CURRENT_SOURCE_DIR}/pythia8_example.cc )
    if (PYTHIA8_VERSION VERSION_GREATER "8.300")
      target_compile_definitions(pythia8_example.exe PRIVATE "-DHEPMC3_USE_INTERFACE_FROM_PYTHIA8")
      target_include_directories(pythia8_example.exe PRIVATE ${PYTHIA8_INCLUDE_DIR})
    else()
      target_include_directories(pythia8_example.exe PRIVATE ${PROJECT_SOURCE_DIR}/interfaces/pythia8/include/Pythia8  ${PYTHIA8_INCLUDE_DIR})
    endif()
   endif()
endif()


if(Pythia8_FOUND)
  target_link_libraries(pythia8_example.exe ${PYTHIA8_LIBRARIES} ${HEPMC3_LIB})
  set_target_properties(pythia8_example.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR})

# create environment scripts

  set(bindir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")
  get_filename_component(PYTHIA8_LIB_DIR ${PYTHIA8_LIBRARY} PATH)
  set(libdir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:${PYTHIA8_LIB_DIR}")
  set(pythia8data "${PYTHIA8_XMLDOC_DIR}")

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pythia8_example_env.sh.in
               ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/pythia8_example_env.sh
               @ONLY
               )
# installs
  if (USE_INSTALLED_HEPMC3)
    install(TARGETS pythia8_example.exe DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(FILES pythia8_ee_to_Z_to_tautau.conf DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/ DESTINATION ${CMAKE_INSTALL_BINDIR}
           FILES_MATCHING PATTERN "*env.sh*")
  endif()
else()
  message(STATUS "HepMC3 examples: Pythia8 package not found. Please check if PYTHIA8_ROOT_DIR is set properly. pythia8_example disabled")
endif()
