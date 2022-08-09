
set(basic_sources   ${CMAKE_CURRENT_SOURCE_DIR}/convert_example.cc
                    ${CMAKE_CURRENT_SOURCE_DIR}/cmdline.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/WriterDOT.cc
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/WriterHEPEVTZEUS.cc )
set(basic_includes ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} )
set(basic_libraries ${HEPMC3_LIB})
set(basic_defines   HEPMCCONVERT_EXTENSION_HEPEVTZEUS=1   HEPMCCONVERT_EXTENSION_DOT=1)
set(basic_dependencies )
set(basic_libdir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" )
find_package(ZLIB)
if(ZLIB_FOUND)
  message(STATUS "HepMC3 ConvertExample: zlib found. zlib compressed files support enabled. ${ZLIB_INCLUDE_DIR} ${ZLIB_LIBRARIES} ")
  list(APPEND basic_includes ${ZLIB_INCLUDE_DIR})
  list(APPEND basic_libraries ${ZLIB_LIBRARIES})
  list(APPEND basic_defines HEPMC3_USE_COMPRESSION)
  list(APPEND basic_defines HEPMC3_Z_SUPPORT=1)
  # create environment scripts
  #set(basic_libdir "${basic_libdir}:${ZLIB_ROOT}")
elseif()
  message(STATUS "HepMC3 ConvertExample: zlib not found. zlib compressed files support disabled.")
endif()
find_package(LibLZMA)
if(LIBLZMA_FOUND)
  message(STATUS "HepMC3 ConvertExample: lzma found. lzma compressed files support enabled. ${LIBLZMA_INCLUDE_DIRS} ${LIBLZMA_LIBRARIES} ")
  list(APPEND basic_includes ${LIBLZMA_INCLUDE_DIRS})
  list(APPEND basic_libraries ${LIBLZMA_LIBRARIES})
  list(APPEND basic_defines HEPMC3_USE_COMPRESSION)
  list(APPEND basic_defines HEPMC3_LZMA_SUPPORT=1)
elseif()
  message(STATUS "HepMC3 ConvertExample: lzma not found. lzma compressed files support disabled.")
endif()
find_package(BZip2)
if(BZIP2_FOUND)
  message(STATUS "HepMC3 ConvertExample: bzip2 found. bzip2 compressed files support enabled. ${BZIP2_INCLUDE_DIRS} ${BZIP2_LIBRARIES} ")
  list(APPEND basic_includes ${BZIP2_INCLUDE_DIRS})
  list(APPEND basic_libraries ${BZIP2_LIBRARIES})
  list(APPEND basic_defines HEPMC3_USE_COMPRESSION)
  list(APPEND basic_defines HEPMC3_BZ2_SUPPORT=1)
elseif()
  message(STATUS "HepMC3 ConvertExample: bzip not found. bzip2 compressed files support disabled.")
endif()
if(ROOT_FOUND)
  list(APPEND basic_includes ${PROJECT_SOURCE_DIR}/rootIO/include  ${ROOT_INCLUDE_DIRS})
  list(APPEND basic_sources ${CMAKE_CURRENT_SOURCE_DIR}/src/WriterRootTreeOPAL.cc )
  list(APPEND basic_libraries ${ROOT_LIBRARIES} ${HEPMC3_ROOTIO_LIB})
  list(APPEND basic_defines HEPMCCONVERT_EXTENSION_ROOTTREEOPAL=1)
  list(APPEND basic_dependencies dictionary)
  # create environment scripts
  set(basic_libdir "${basic_libdir}:${ROOT_LIBRARY_DIR}")
endif()
if(Protobuf_FOUND)
  list(APPEND basic_includes ${PROJECT_SOURCE_DIR}/protobufIO/include  ${Protobuf_INCLUDE_DIRS})
  list(APPEND basic_libraries ${Protobuf_LIBRARIES} ${HEPMC3_PROTOBUFIO_LIB})
  # create environment scripts
  set(basic_libdir "${basic_libdir}:${Protobuf_LIBRARY_DIR}")
endif()


if(${CMAKE_VERSION} VERSION_GREATER "3.8")
########################################################################
#FindPython is not reenterant, so we have to hack the cache.
########################################################################
  unset(_Python3_LIBRARY_RELEASE CACHE) # Attempt to make FindPython reenterant
  unset(_Python3_NumPy_INCLUDE_DIR CACHE) # Attempt to make FindPython reenterant
  find_package(Python3 3.5 COMPONENTS Development.Embed NumPy)
  if (Python3_FOUND)
     list(APPEND basic_includes  ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
     list(APPEND basic_sources   ${CMAKE_CURRENT_SOURCE_DIR}/src/ReaderuprootTree.cc )
     list(APPEND basic_libraries ${Python3_LIBRARIES})
     list(APPEND basic_defines HEPMCCONVERT_EXTENSION_UPROOTTREEREADER=1)
  elseif()
     message(STATUS "HepMC3 ConvertExample: Python3.5 with NumPy not found. uproot support disabled.")
  endif()
  unset(_Python3_LIBRARY_RELEASE CACHE) # Attempt to make FindPython reenterant
  unset(_Python3_NumPy_INCLUDE_DIR CACHE) # Attempt to make FindPython reenterant
########################################################################
endif()

add_executable(convert_example.exe  ${basic_sources} )
add_library(Analysis SHARED src/AnalysisExample.cc)
target_link_libraries(Analysis ${HEPMC3_LIB} )

if (MSVC)
  set(  MSVCLINKFLAGS "/WHOLEARCHIVE:${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/HepMC3-static.lib")
  if(ZLIB_FOUND)
    set(  MSVCLINKFLAGS "${MSVCLINKFLAGS} /WHOLEARCHIVE:${ZLIB_LIBRARIES}")
  endif()
  set_target_properties(convert_example.exe PROPERTIES LINK_FLAGS ${MSVCLINKFLAGS})
endif()
include_directories(${basic_includes})
target_link_libraries(convert_example.exe ${basic_libraries} )
target_compile_definitions(convert_example.exe PRIVATE ${basic_defines} )
if (basic_dependencies)
  add_dependencies( convert_example.exe ${basic_dependencies} )
endif()
# create environment scripts
set(bindir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")
set(libdir ${basic_libdir})
set_target_properties(convert_example.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/$<0:>)
set_target_properties(Analysis PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/$<0:>)

configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/convert_example_env.sh.in
	    ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/convert_example_env.sh @ONLY)

# installs
if (USE_INSTALLED_HEPMC3)
  install(TARGETS convert_example.exe DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(TARGETS Analysis DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${CMAKE_INSTALL_BINDIR} FILES_MATCHING PATTERN "*.pcm" PATTERN "*config"  PATTERN "CMakeFiles" PATTERN  "*.ilk"  EXCLUDE)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/ DESTINATION ${CMAKE_INSTALL_BINDIR}  FILES_MATCHING PATTERN "*env.sh*" )
endif()

#This is an option for developers only
if ( HEPMC3_USE_GENGETOPT )
  find_program(GENGETOPT gengetopt)
  if( NOT GENGETOPT )
    message( STATUS  "HepMC3 examples: gengetopt not found. Old scripts for options parsing will be used." )
  else()
    message( STATUS  "HepMC3 examples: gengetopt found in ${GENGETOPT}. New scripts for options parsing can be generated." )
    add_custom_command(
      OUTPUT  ${CMAKE_CURRENT_SOURCE_DIR}/cmdline.c ${CMAKE_CURRENT_SOURCE_DIR}/cmdline.h
      COMMAND ${GENGETOPT} -N --no-version --input=${CMAKE_CURRENT_SOURCE_DIR}/cmdline.ggo
      --file-name cmdline --output-dir=${CMAKE_CURRENT_SOURCE_DIR}
      --include-getopt --show-required --unnamed-opts=\"inputfile outputfile\"
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cmdline.ggo
      COMMENT        "HepMC3 examples: Generating getopt parser code..."
    )
  endif()
  add_custom_target( cmdlineoptions DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cmdline.c ${CMAKE_CURRENT_SOURCE_DIR}/cmdline.h)
  add_dependencies(  convert_example.exe cmdlineoptions )
endif()


