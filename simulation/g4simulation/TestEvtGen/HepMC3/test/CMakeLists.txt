macro( hepmc_test testname )
  message( STATUS "HepMC3 test: building ${testname} for ${testname}" )
  find_file( ${testname}_source ${testname}.cc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
  ADD_EXECUTABLE(${testname} ${${testname}_source}   ${ARGN} )
  target_include_directories(${testname} BEFORE
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src   ${CMAKE_CURRENT_SOURCE_DIR}/../include/)
  target_link_libraries(${testname} PUBLIC  PRIVATE HepMC3 )
  set(srcdir ${CMAKE_CURRENT_SOURCE_DIR} )
  set(builddir ${CMAKE_CURRENT_BINARY_DIR} )
  add_test( ${testname} ${CMAKE_CURRENT_BINARY_DIR}/${testname} )
  SET_TESTS_PROPERTIES( ${testname} PROPERTIES ENVIRONMENT "PATH=${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/;$ENV{PATH}")
  file(GLOB INPUT${testname}  "${CMAKE_CURRENT_SOURCE_DIR}/input${testname}*" )
  file(COPY ${INPUT${testname}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  string(REPLACE "test" "" ${testname} name)
  file(GLOB INPUT${name}  "${CMAKE_CURRENT_SOURCE_DIR}/input${name}*" )
  file(COPY ${INPUT${name}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endmacro( hepmc_test )

macro( hepmc_root_test testname )
  message( STATUS "HepMC3 test: building ROOT test ${testname} " )
  find_file( ${testname}_source ${testname}.cc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
  ADD_EXECUTABLE(${testname} ${${testname}_source}  ${ARGN} )
  target_include_directories(${testname} BEFORE
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src   ${CMAKE_CURRENT_SOURCE_DIR}/../include/ ${CMAKE_CURRENT_SOURCE_DIR}/../rootIO/include/ )
  target_link_libraries(${testname} PRIVATE HepMC3 HepMC3rootIO)

  set(srcdir ${CMAKE_CURRENT_SOURCE_DIR} )
  set(builddir ${CMAKE_CURRENT_BINARY_DIR} )
  add_test( ${testname} ${testname} )
  SET_TESTS_PROPERTIES( ${testname} PROPERTIES ENVIRONMENT "PATH=${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/;$ENV{PATH}"
  ENVIRONMENT "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/;$ENV{LD_LIBRARY_PATH}"
  ENVIRONMENT "DYLD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/;$ENV{DYLD_LIBRARY_PATH}"
  )
  file(GLOB INPUT${testname}  "${CMAKE_CURRENT_SOURCE_DIR}/input${testname}*" )
  file(COPY ${INPUT${testname}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  string(REPLACE "test" "" ${testname} name)
  file(GLOB INPUT${name}  "${CMAKE_CURRENT_SOURCE_DIR}/input${name}*" )
  file(COPY ${INPUT${name}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endmacro( hepmc_root_test )


macro( hepmc_protobuf_test testname )
  message( STATUS "HepMC3 test: building protobuf test ${testname} " )
  find_file( ${testname}_source ${testname}.cc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
  ADD_EXECUTABLE(${testname} ${${testname}_source}  ${ARGN} )
  target_include_directories(${testname} BEFORE
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src   ${CMAKE_CURRENT_SOURCE_DIR}/../include/ ${CMAKE_CURRENT_SOURCE_DIR}/../protobufIO/include/ )
  target_link_libraries(${testname} PRIVATE HepMC3 HepMC3protobufIO)

  set(srcdir ${CMAKE_CURRENT_SOURCE_DIR} )
  set(builddir ${CMAKE_CURRENT_BINARY_DIR} )
  add_test( ${testname} ${testname} )
  SET_TESTS_PROPERTIES( ${testname} PROPERTIES ENVIRONMENT "PATH=${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/;$ENV{PATH}")
  file(GLOB INPUT${testname}  "${CMAKE_CURRENT_SOURCE_DIR}/input${testname}*" )
  file(COPY ${INPUT${testname}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  string(REPLACE "test" "" ${testname} name)
  file(GLOB INPUT${name}  "${CMAKE_CURRENT_SOURCE_DIR}/input${name}*" )
  file(COPY ${INPUT${name}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endmacro( hepmc_protobuf_test )


set( CMAKE_INCLUDE_DIRECTORIES_BEFORE ON )
cmake_policy(SET CMP0057 NEW)

option(HEPMC3_TEST_MCTESTER  "Use MC-TESTER in tests if available" ON)
option(HEPMC3_TEST_PHOTOSPP  "Use Photospp in tests if available" ON)
option(HEPMC3_TEST_TAUOLAPP  "Use Tauola in tests if available" ON)
option(HEPMC3_TEST_PYTHIA8   "Use Pythia8 in tests if available" ON)
option(HEPMC3_TEST_THREADS   "Use Threads in tests if available" ON)
option(HEPMC3_TEST_HEPMC2    "Use HepMC2 in tests if available" ON)
option(HEPMC3_TEST_VALGRIND  "Use valgrind in tests if available" OFF)
option(HEPMC3_TEST_ZLIB      "Use zlib in tests if available" ON)
option(HEPMC3_TEST_LZMA      "Use lzma in tests if available" ON)
option(HEPMC3_TEST_BZIP2     "Use bzip2 in tests if available" ON)

set(HEPMC3_TEST_PACKAGES_LIST "")
if (HEPMC3_TEST_MCTESTER)
 LIST(APPEND HEPMC3_TEST_PACKAGES_LIST "MCTester")
endif()
if (HEPMC3_TEST_PHOTOSPP)
 LIST(APPEND HEPMC3_TEST_PACKAGES_LIST "Photospp")
endif()
if (HEPMC3_TEST_TAUOLAPP)
 LIST(APPEND HEPMC3_TEST_PACKAGES_LIST "Tauolapp")
endif()
if (HEPMC3_TEST_PYTHIA8)
 LIST(APPEND HEPMC3_TEST_PACKAGES_LIST "Pythia8")
endif()
if (HEPMC3_TEST_THREADS)
 LIST(APPEND HEPMC3_TEST_PACKAGES_LIST "Threads")
endif()
if (HEPMC3_TEST_HEPMC2)
 LIST(APPEND HEPMC3_TEST_PACKAGES_LIST "HepMC2")
endif()
if (HEPMC3_TEST_VALGRIND)
 LIST(APPEND HEPMC3_TEST_PACKAGES_LIST "valgrind")
endif()
if (HEPMC3_TEST_ZLIB)
 LIST(APPEND HEPMC3_TEST_PACKAGES_LIST "zlib")
endif()
if (HEPMC3_TEST_LZMA)
 LIST(APPEND HEPMC3_TEST_PACKAGES_LIST "lzma")
endif()
if (HEPMC3_TEST_BZIP2)
 LIST(APPEND HEPMC3_TEST_PACKAGES_LIST "bzip2")
endif()
message(STATUS "HepMC3 test: The following packages are enabled for tests HEPMC3_TEST_PACKAGES_LIST=${HEPMC3_TEST_PACKAGES_LIST}" )

option(HEPMC3_USE_INTERFACE_FROM_MCTESTER "Use interface from installed MC-TESTER" ON)
option(HEPMC3_USE_INTERFACE_FROM_PHOTOSPP "Use interface from installed Photos" ON)
option(HEPMC3_USE_INTERFACE_FROM_TAUOLAPP "Use interface from installed Tauola" ON)
option(HEPMC3_USE_INTERFACE_FROM_PYTHIA8  "Use interface from installed Pythia8" ON)

file(GLOB INPUTFILES  "${CMAKE_CURRENT_SOURCE_DIR}/*.input" )
file(COPY ${INPUTFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set( HepMC_tests
        testBoost
        testFourVector
        testMass
        testPrintBug
        testPolarization
        testUnits
        testMultipleCopies
        testWeights
        testLoops
        testDelete
        testSkip1
        testIO1
        testIO3
        testIO6
        testIO7
        testIO8
        testReaderFactory1
        testSingleVertexHepMC2
        testAttributes
        testHEPEVTWrapper1

         )
set( HepMC_root_tests
        testIO2
        testIO4
        testReaderFactory2
        testRoot300
        testRootTree300
        )

set( HepMC_search_tests

        )
if (WIN32)
  set( HepMC_protobuf_tests )
else()
  set( HepMC_protobuf_tests
        testIO20
        testIO21
        testIO22
        testIO23
        testIO24
        testIO25
        testIO26
        testIO27
        testIO28
        testIO29
        )
endif()

set(srcdir ${CMAKE_CURRENT_SOURCE_DIR} )
if ("Threads" IN_LIST HEPMC3_TEST_PACKAGES_LIST)
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads)
endif()
if("Threads" IN_LIST HEPMC3_TEST_PACKAGES_LIST AND Threads_FOUND AND NOT (${CMAKE_CXX_COMPILER_ID} MATCHES "PGI") AND NOT (${CMAKE_CXX_COMPILER_ID} MATCHES "NVHPC"))
  message(STATUS "HepMC3 test: Threads library found.  Thread safety tests enabled")
  list( APPEND HepMC_tests "testThreads1" )
  list( APPEND HepMC_tests "testIO10" )
  if (HEPMC3_ENABLE_SEARCH)
    list( APPEND HepMC_search_tests "testThreadssearch" )
  endif()
else()
  message(STATUS "HepMC3 test: No threads library found or disabled for compiller with ID ${CMAKE_CXX_COMPILER_ID}. Thread safety tests are disabled")
endif()

set(compress_tests "testIO9")
if ("zlib" IN_LIST HEPMC3_TEST_PACKAGES_LIST)
  find_package(ZLIB)
  if(ZLIB_FOUND)
    message(STATUS "HepMC3 test: ZLIB library found. zlib tests enabled")
    list( APPEND HepMC_tests ${compress_tests} )
  endif()
endif()
if ("lzma" IN_LIST HEPMC3_TEST_PACKAGES_LIST)
  find_package(LibLZMA)
  if(LIBLZMA_FOUND)
    message(STATUS "HepMC3 test: LZMA library found.  lzma tests enabled")
    list( APPEND HepMC_tests ${compress_tests} )
  endif()
endif()
if ("bzip2" IN_LIST HEPMC3_TEST_PACKAGES_LIST)
  find_package(BZip2)
  if(LIBLZMA_FOUND)
    message(STATUS "HepMC3 test: BZip2 library found.  BZip2 tests enabled")
    list( APPEND HepMC_tests ${compress_tests} )
  endif()
endif()  
## Here we add the tests
list(REMOVE_DUPLICATES HepMC_tests)
foreach ( test ${HepMC_tests} )
  hepmc_test( ${test} )
  set_tests_properties(${test} PROPERTIES TIMEOUT 2000)  #Relevant for testLoops
endforeach ( test ${HepMC_tests} )

if(HEPMC3_ENABLE_SEARCH )
  foreach ( test ${HepMC_search_tests} )
    hepmc_test( ${test} )
  target_link_libraries(${test} PRIVATE HepMC3search)
  target_include_directories(${test} BEFORE PRIVATE  ${CMAKE_CURRENT_SOURCE_DIR}/../search/include)
  endforeach ( test ${HepMC_search_tests} )
endif()

if(HEPMC3_ENABLE_PROTOBUFIO )
  foreach ( test ${HepMC_protobuf_tests} )
    hepmc_protobuf_test( ${test} )
  target_link_libraries(${test} PRIVATE HepMC3protobufIO)
  target_include_directories(${test} BEFORE PRIVATE  ${CMAKE_CURRENT_SOURCE_DIR}/../protobufIO/include)
  endforeach ( test ${HepMC_protobuf_tests} )
endif()


## Treatment of special needs in the HepMC tests
if(Threads_FOUND AND NOT (${CMAKE_CXX_COMPILER_ID} MATCHES "PGI") AND NOT (${CMAKE_CXX_COMPILER_ID} MATCHES "NVHPC") )
  if(THREADS_HAVE_PTHREAD_ARG)
    target_compile_options(testThreads1 PUBLIC "-pthread")
    target_compile_options(testIO10 PUBLIC "-pthread")
    if (HEPMC3_ENABLE_SEARCH)
     target_compile_options(testThreadssearch PUBLIC "-pthread")
    endif()    
  endif()
  if(CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(testThreads1 PUBLIC "${CMAKE_THREAD_LIBS_INIT}")
    target_link_libraries(testIO10 PUBLIC "${CMAKE_THREAD_LIBS_INIT}")
    if (HEPMC3_ENABLE_SEARCH)
     target_link_libraries(testThreadssearch PUBLIC "${CMAKE_THREAD_LIBS_INIT}")
    endif() 
  endif()
endif()
if (WIN32)
  if ("zlib" IN_LIST HEPMC3_TEST_PACKAGES_LIST AND ZLIB_FOUND)
     get_filename_component(z ${ZLIB_LIBRARIES} DIRECTORY)
     find_file(zdll NAME zlib.dll PATHS ${z} ${z}/../bin/)
     if (zdll)
       file(COPY ${zdll} DESTINATION ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/)
     endif()
   endif()
   if ("lzma" IN_LIST HEPMC3_TEST_PACKAGES_LIST AND LIBLZMA_FOUND)
     get_filename_component(l ${LIBLZMA_LIBRARIES} DIRECTORY)
     find_file(lzmadll NAME lzma.dll PATHS ${l} ${l}/../bin/)
     if (lzmadll)
       file(COPY ${lzmadll} DESTINATION ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/)
     endif()
   endif()
   if ("bzip2" IN_LIST HEPMC3_TEST_PACKAGES_LIST AND BZIP2_FOUND)
     get_filename_component(b ${BZIP2_LIBRARIES} DIRECTORY)
     find_file(bz2dll NAME libbz2.dll PATHS ${b} ${b}/../bin/)
     if (bz2dll)
       file(COPY ${bz2dll} DESTINATION ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/)
     endif()
   endif()
endif()

foreach ( ctest ${compress_tests} )
  if ("zlib" IN_LIST HEPMC3_TEST_PACKAGES_LIST AND ZLIB_FOUND)
    target_compile_options(${ctest} PUBLIC "-DHEPMC3_USE_COMPRESSION")
    target_compile_options(${ctest} PUBLIC "-DHEPMC3_Z_SUPPORT=1")
    target_link_libraries(${ctest} PRIVATE "${ZLIB_LIBRARIES}" )
    target_include_directories(${ctest} PRIVATE "${ZLIB_INCLUDE_DIRS}" )
  endif()
  if ("lzma" IN_LIST HEPMC3_TEST_PACKAGES_LIST AND LIBLZMA_FOUND)
    target_compile_options(${ctest} PUBLIC "-DHEPMC3_USE_COMPRESSION")
    target_compile_options(${ctest} PUBLIC "-DHEPMC3_LZMA_SUPPORT=1")
    target_link_libraries(${ctest} PRIVATE "${LIBLZMA_LIBRARIES}" )
    target_include_directories(${ctest} PRIVATE "${LIBLZMA_INCLUDE_DIRS}" )
  endif()
  if ("bzip2" IN_LIST HEPMC3_TEST_PACKAGES_LIST AND BZIP2_FOUND)
    target_compile_options(${ctest} PUBLIC "-DHEPMC3_USE_COMPRESSION")
    target_compile_options(${ctest} PUBLIC "-DHEPMC3_BZ2_SUPPORT=1")
    target_link_libraries(${ctest} PRIVATE "${BZIP2_LIBRARIES}" )
    target_include_directories(${ctest} PRIVATE "${BZIP2_INCLUDE_DIRS}" )
  endif()
endforeach ( ctest ${compress_tests} )
#So far the test with ROOT have issues on WIN32
if(HEPMC3_ENABLE_ROOTIO AND (NOT WIN32) )
  foreach ( test ${HepMC_root_tests} )
    hepmc_root_test( ${test} )
  endforeach ( test ${HepMC_root_tests} )
  target_link_libraries(testIO4 PRIVATE ${ROOT_Hist_LIBRARY})
endif()


set( libMCEG_sources ValidationControl.cc   SimpleEventTool.cc OutputValidationTool.cc)
set( libMCEG_sources_HEPMC3  )
set( libMCEG_link_libraries_HEPMC2  )
set( libMCEG_link_libraries_HEPMC3  )
set( libMCEG_include_directories   ${CMAKE_CURRENT_SOURCE_DIR}/../include/  )
set( MCEG_tests )

set( convert_tests )


if ("Pythia8" IN_LIST HEPMC3_TEST_PACKAGES_LIST)
    find_package(Pythia8 8.220)
endif()
if(NOT Pythia8_FOUND)
  message(STATUS "HepMC3 test: Pythia8 package not found, was disabled or too old. The minimal supported Pythia8 version is 8.220 (8.300 if the Pythia8 interface is used). Please check if PYTHIA8_ROOT_DIR is set properly.")
else()
  message(STATUS "HepMC3 test: Pythia8 package found. ${PYTHIA8_LIBRARIES} ${PYTHIA8_VERSION}")
  list( APPEND libMCEG_sources    "PythiaValidationTool.cc" )
  list( APPEND libMCEG_link_libraries_HEPMC2  "${PYTHIA8_LIBRARIES}")
  list( APPEND libMCEG_link_libraries_HEPMC3  "${PYTHIA8_LIBRARIES}")
  if (PYTHIA8_VERSION VERSION_GREATER "8.300")
    list( APPEND libMCEG_compile_definitions  "-DHEPMC3_USE_INTERFACE_FROM_PYTHIA8")
  else()
    message(STATUS "HepMC3 test: WARNING: Will use deprecated built-in HepMC3 interface for Pythia8.")
    list( APPEND libMCEG_include_directories  "${PROJECT_SOURCE_DIR}/interfaces/pythia8/include/Pythia8")
  endif()
  list( APPEND libMCEG_include_directories  "${PYTHIA8_INCLUDE_DIRS}")
  list( APPEND libMCEG_compile_definitions  "-DPYTHIA8")

  list( APPEND MCEG_tests "testPythia1" )
  list( APPEND MCEG_tests "testPythia2" )
  list( APPEND MCEG_tests "testPythia3" )
endif()

if ("Photospp" IN_LIST HEPMC3_TEST_PACKAGES_LIST  AND Pythia8_FOUND)
   find_package(Photospp COMPONENTS HEPMC2 HEPMC3)
endif()
if(NOT Photospp_FOUND)
  message(STATUS "HepMC3 test: Photospp package not found or disabled. Please check if PHOTOSPP_ROOT_DIR is set properly.")
else()
  message(STATUS "HepMC3 test: Photospp package found. Building Photospp tests. PHOTOSPP_VERSION ${PHOTOSPP_VERSION} PHOTOSPP_LIBRARIES ${PHOTOSPP_LIBRARIES} PHOTOSPP_HEPMC2_LIBRARIES ${PHOTOSPP_HEPMC2_LIBRARIES} PHOTOSPP_HEPMC3_LIBRARIES ${PHOTOSPP_HEPMC3_LIBRARIES} PHOTOSPP_INCLUDE_DIRS ${PHOTOSPP_INCLUDE_DIRS}")
  list( APPEND  libMCEG_sources   "PhotosValidationTool.cc" )
  list( APPEND  libMCEG_link_libraries_HEPMC2  "${PHOTOSPP_LIBRARIES}"  "${PHOTOSPP_HEPMC2_LIBRARIES}")
  list( APPEND  libMCEG_link_libraries_HEPMC3  "${PHOTOSPP_LIBRARIES}" )
  list( APPEND  libMCEG_link_libraries_HEPMC3  "${PHOTOSPP_HEPMC3_LIBRARIES}" )
  list( APPEND  libMCEG_include_directories  "${PHOTOSPP_INCLUDE_DIRS}")
  list( APPEND  libMCEG_compile_definitions  "-DPHOTOSPP")
  list( APPEND  MCEG_tests "testPhotos1" )
endif()

if ("Tauolapp" IN_LIST HEPMC3_TEST_PACKAGES_LIST AND Pythia8_FOUND)
  find_package(Tauolapp COMPONENTS HEPMC2 HEPMC3)
endif()
if(NOT Tauolapp_FOUND)
  message(STATUS "HepMC3 test: Tauolapp package not found or disabled. Please check if TAUOLA_ROOT_DIR is set properly.")
else()
  message(STATUS "HepMC3 test: Tauolapp package found. Building Tauolapp tests. TAUOLAPP_VERSION ${TAUOLAPP_VERSION} TAUOLAPP_LIBRARIES ${TAUOLAPP_LIBRARIES} TAUOLAPP_HEPMC2_LIBRARIES ${TAUOLAPP_HEPMC2_LIBRARIES} TAUOLAPP_HEPMC3_LIBRARIES ${TAUOLAPP_HEPMC3_LIBRARIES} TAUOLAPP_INCLUDE_DIRS ${TAUOLAPP_INCLUDE_DIRS}")
  list( APPEND  libMCEG_sources "TauolaValidationTool.cc" )
  list( APPEND  libMCEG_link_libraries_HEPMC2  "${TAUOLAPP_LIBRARIES}" "${TAUOLAPP_HEPMC2_LIBRARIES}" )
  list( APPEND  libMCEG_link_libraries_HEPMC3  "${TAUOLAPP_LIBRARIES}" )
  list( APPEND  libMCEG_link_libraries_HEPMC3  "${TAUOLAPP_HEPMC3_LIBRARIES}" )
  list( APPEND  libMCEG_include_directories  "${TAUOLAPP_INCLUDE_DIRS}")
  list( APPEND  libMCEG_compile_definitions  "-DTAUOLAPP")
  list( APPEND   MCEG_tests "testTauola1" )
endif()

if(HEPMC3_ENABLE_ROOTIO)
  if ("MCTester" IN_LIST HEPMC3_TEST_PACKAGES_LIST)
    if (HEPMC3_USE_INTERFACE_FROM_MCTESTER)
      find_package(MCTester COMPONENTS HEPMC2 HEPMC3)
    else()
      find_package(MCTester COMPONENTS HEPMC2)
    endif()
  endif()
  if(NOT MCTester_FOUND)
    message(STATUS "HepMC3 test: MCTester package not found or disabled. Please check if MCTESTER_ROOT_DIR is set properly.")
  else()
    message(STATUS "HepMC3 test: MCTester package found. Building MCTester tests. MCTESTER_VERSION ${MCTESTER_VERSION} MCTESTER_LIBRARIES ${MCTESTER_LIBRARIES} MCTESTER_HEPMC2_LIBRARIES ${MCTESTER_HEPMC2_LIBRARIES} MCTESTER_HEPMC3_LIBRARIES ${MCTESTER_HEPMC3_LIBRARIES} MCTESTER_INCLUDE_DIRS ${MCTESTER_INCLUDE_DIRS}")
    list( APPEND  libMCEG_sources "McTesterValidationTool.cc" )
    list( APPEND  libMCEG_link_libraries_HEPMC2  "${MCTESTER_LIBRARIES}" "${ROOT_LIBRARIES}" "${MCTESTER_HEPMC2_LIBRARIES}" )
    list( APPEND  libMCEG_link_libraries_HEPMC3  "${MCTESTER_LIBRARIES}" "${ROOT_LIBRARIES}" )
    list( APPEND  libMCEG_link_libraries_HEPMC3  "${MCTESTER_HEPMC3_LIBRARIES}")
    list( APPEND  libMCEG_include_directories    "${MCTESTER_INCLUDE_DIRS}" "${ROOT_INCLUDE_DIRS}")
    list( APPEND  libMCEG_compile_definitions  "-DMCTESTER")
    list( APPEND  MCEG_tests "testMCTester1" )
  endif()
endif()


add_library(MCEG3 SHARED ${libMCEG_sources} ${libMCEG_sources_HEPMC3})
target_link_libraries(MCEG3  PRIVATE HepMC3 PUBLIC ${libMCEG_link_libraries_HEPMC3} )
target_compile_definitions(MCEG3  PRIVATE ${libMCEG_compile_definitions})
target_include_directories(MCEG3   BEFORE
PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
PRIVATE ${libMCEG_include_directories}
PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src   ${CMAKE_CURRENT_SOURCE_DIR}/../include/)



macro( hepmc_mceghepmc3_test testname )
  message( STATUS "HepMC3 test: building MCEG3 test ${testname} with HEPMC3" )
  find_file( ${testname}_source ${testname}.cc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
  ADD_EXECUTABLE(${testname} ${${testname}_source}  ${ARGN} )
  target_include_directories(${testname} BEFORE
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src   ${CMAKE_CURRENT_SOURCE_DIR}/../include/   )
  #The flag is needed for Taoula that conflates generic C++ interface with the HepMC2 interface
  #FIXME put here Visual studio flags when needed
  #For Darwin target_link_libraries(${testname}  "-undefined dynamic_lookup"  HepMC3  MCEG3)
if (APPLE)
  target_link_libraries(${testname}  "-undefined dynamic_lookup"  HepMC3  MCEG3)
else()
  target_link_libraries(${testname}  "-Wl,--allow-shlib-undefined"  HepMC3  MCEG3)
endif()
  set(srcdir ${CMAKE_CURRENT_SOURCE_DIR} )
  set(builddir ${CMAKE_CURRENT_BINARY_DIR} )
  add_test( ${testname} ${testname} )
  SET_TESTS_PROPERTIES( ${testname} PROPERTIES ENVIRONMENT "PATH=${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/;$ENV{PATH}" ENVIRONMENT "PYTHIA8DATA=${PYTHIA8_XMLDOC_DIR}")
  file(GLOB INPUT${testname}  "${CMAKE_CURRENT_SOURCE_DIR}/input${testname}*" )
  file(COPY ${INPUT${testname}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  string(REPLACE "test" "" ${testname} name)
  file(GLOB INPUT${name}  "${CMAKE_CURRENT_SOURCE_DIR}/input${name}*" )
  file(COPY ${INPUT${name}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endmacro( hepmc_mceghepmc3_test )


foreach ( test ${MCEG_tests} )
  hepmc_mceghepmc3_test( ${test} )
endforeach ( test ${MCEG_tests} )




if ("HepMC2" IN_LIST HEPMC3_TEST_PACKAGES_LIST)
  find_package(HepMC2 QUIET)
endif()

if(NOT HepMC2_FOUND)
  message(STATUS "HepMC3 test: HEPMC2 package not found or disabled. Please check if HEPMC2_ROOT_DIR is set properly.")
else()
  message(STATUS "HepMC3 test: HEPMC2 package found. HEPMC2_VERSION ${HEPMC2_VERSION} HEPMC2_LIBRARIES ${HEPMC2_LIBRARIES} HEPMC2_INCLUDE_DIRS ${HEPMC2_INCLUDE_DIRS}")
  add_library(MCEG2 SHARED ${libMCEG_sources})
  target_link_libraries(MCEG2  PRIVATE ${HEPMC2_LIBRARIES} PUBLIC ${libMCEG_link_libraries_HEPMC2} )
  target_compile_definitions(MCEG2  PRIVATE ${libMCEG_compile_definitions} -DHEPMC2)
  target_include_directories(MCEG2   BEFORE
  PUBLIC ${HEPMC2_INCLUDE_DIRS}
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${libMCEG_include_directories}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src   )

  list( APPEND convert_tests "testConvert1" )
  list( APPEND convert_tests "testConvert2" )

macro( hepmc_mceghepmc2_test testname )
  message( STATUS "HepMC3 test: building MCEG2 test ${testname} with HEPMC2" )
  find_file( ${testname}_source ${testname}.cc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
  ADD_EXECUTABLE(hepmc2_${testname} ${${testname}_source}  ${ARGN} )
  target_compile_definitions(hepmc2_${testname}  PRIVATE  -DHEPMC2)
  target_include_directories(${testname} BEFORE
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src   ${CMAKE_CURRENT_SOURCE_DIR}/../include/   )
  target_link_libraries(hepmc2_${testname} PRIVATE ${HEPMC2_LIBRARIES}   MCEG2   )
  set(srcdir ${CMAKE_CURRENT_SOURCE_DIR} )
  set(builddir ${CMAKE_CURRENT_BINARY_DIR} )
  add_test( hepmc2_${testname} hepmc2_${testname} )
  SET_TESTS_PROPERTIES(  hepmc2_${testname} PROPERTIES ENVIRONMENT "PATH=${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/;$ENV{PATH}" ENVIRONMENT "PYTHIA8DATA=${PYTHIA8_XMLDOC_DIR}")
  file(GLOB INPUThepmc2_${testname}  "${CMAKE_CURRENT_SOURCE_DIR}/input${testname}*" )
  file(COPY ${INPUThepmc2_${testname}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  string(REPLACE "test" "" ${testname} name)
  file(GLOB INPUThepmc2_${name}  "${CMAKE_CURRENT_SOURCE_DIR}/input${name}*" )
  file(COPY ${INPUThepmc2_${name}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endmacro( hepmc_mceghepmc2_test )



foreach ( test ${MCEG_tests} )
  hepmc_mceghepmc2_test( ${test} )
endforeach ( test ${MCEG_tests} )

macro( hepmc_convert_test testname )
  message( STATUS "HepMC3 test: building Convert test ${testname} with HEPMC2" )
  find_file( ${testname}_source ${testname}.cc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
  message( STATUS "HepMC3 test: building Convert test ${testname} ${${testname}_source} with HEPMC2" )
  ADD_EXECUTABLE(convert_${testname} ${${testname}_source}  ${ARGN} )
  target_include_directories(convert_${testname} BEFORE
  PUBLIC ${HEPMC2_INCLUDE_DIRS}
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src   ${CMAKE_CURRENT_SOURCE_DIR}/../include/   "${PROJECT_SOURCE_DIR}/interfaces/HepMCCompatibility/include" )
  target_link_libraries(convert_${testname} PRIVATE ${HEPMC2_LIBRARIES}   HepMC3  )
  set(srcdir ${CMAKE_CURRENT_SOURCE_DIR} )
  set(builddir ${CMAKE_CURRENT_BINARY_DIR} )
  add_test( convert_${testname} convert_${testname} )
  SET_TESTS_PROPERTIES(  convert_${testname} PROPERTIES ENVIRONMENT "PATH=${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/;$ENV{PATH}")
  file(GLOB INPUTconvert_${testname}  "${CMAKE_CURRENT_SOURCE_DIR}/input${testname}*" )
  file(COPY ${INPUTconvert_${testname}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  string(REPLACE "test" "" ${testname} name)
  file(GLOB INPUTconvert_${name}  "${CMAKE_CURRENT_SOURCE_DIR}/input${name}*" )
  file(COPY ${INPUTconvert_${name}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endmacro( hepmc_convert_test )

foreach ( test ${convert_tests} )
  hepmc_convert_test( ${test} )
endforeach ( test ${convert_tests} )

endif()

#One should use something else for Windows/BSD/MAC


if ("valgrind" IN_LIST HEPMC3_TEST_PACKAGES_LIST)
  find_program( MEMORYCHECK_COMMAND valgrind )
endif()
if(MEMORYCHECK_COMMAND)
  message(STATUS "HepMC3 test: valgrind package found in ${MEMORYCHECK_COMMAND}.  memory check tests enabled.")
  set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full --error-exitcode=1 --default-suppressions=yes --gen-suppressions=all" )
  if (${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
     file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/intelHepMC3.supp" "{\n<HepMC3_at_Linux_with_Intel_1>\nMemcheck:Addr32\nfun:__intel_*\n}\n{\n<HepMC3_at_Linux_with_Intel_2>\nMemcheck:Cond\nfun:*strchr\n}\n"  )
     set( MEMORYCHECK_COMMAND_OPTIONS  "${MEMORYCHECK_COMMAND_OPTIONS} --suppressions=${CMAKE_CURRENT_BINARY_DIR}/IntelHepMC3.supp")
  endif()
elseif()
  message(STATUS "HepMC3 test: valgrind package not found or disabled.  memory check tests disabled.")
endif()


macro(add_memcheck_test binaryname )
  set(memcheck_command "${MEMORYCHECK_COMMAND} ${MEMORYCHECK_COMMAND_OPTIONS}")
  separate_arguments(memcheck_command)
  add_test(memcheck_${binaryname} ${memcheck_command} ./${binaryname} ${ARGN})
  SET_TESTS_PROPERTIES( memcheck_${binaryname} PROPERTIES ENVIRONMENT "PATH=${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/;$ENV{PATH}")
  file(GLOB INPUTmemcheck_${binaryname}  "${CMAKE_CURRENT_SOURCE_DIR}/input${binaryname}*" )
  file(COPY ${INPUTmemcheck_${binaryname}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  string(REPLACE "test" "" ${binaryname} name)
  file(GLOB INPUTmemcheck_${name}  "${CMAKE_CURRENT_SOURCE_DIR}/input${name}*" )
  file(COPY ${INPUTmemcheck_${name}} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endmacro(add_memcheck_test)


if(MEMORYCHECK_COMMAND)
  set( HepMC_memcheck_tests ${HepMC_tests} )
  if(HEPMC3_ENABLE_PROTOBUFIO )
    set( HepMC_memcheck_tests ${HepMC_memcheck_tests} ${HepMC_protobuf_tests} )
  endif()
  list(REMOVE_ITEM HepMC_memcheck_tests testFourVector )
#The story behind is complicated. It seems that on some platforms valgind dislikes the comparison of infinities.
  foreach ( test ${HepMC_memcheck_tests} )
    add_memcheck_test( ${test} )
    set_tests_properties(memcheck_${test} PROPERTIES DEPENDS ${test})
  endforeach ( test ${HepMC_memcheck_tests} )
endif()

