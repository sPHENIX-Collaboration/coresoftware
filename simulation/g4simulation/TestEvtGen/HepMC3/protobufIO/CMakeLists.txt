include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${Protobuf_INCLUDE_DIRS})

protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS HepMC3.proto)

set(HepMC3protobufIO_sources
  src/Writerprotobuf.cc 
  src/Readerprotobuf.cc 
  ${PROTO_SRCS})

include_directories(${CMAKE_CURRENT_BINARY_DIR} 
                    ${CMAKE_CURRENT_BINARY_DIR}/include 
                    ${CMAKE_CURRENT_BINARY_DIR}/include/HepMC3)

#---Define dynamic library---------------------------------------------------
add_library(HepMC3protobufIO SHARED ${HepMC3protobufIO_sources})

target_link_libraries(HepMC3protobufIO HepMC3 protobuf::libprotobuf)

set_target_properties(HepMC3protobufIO 
  PROPERTIES
    SOVERSION 3
    POSITION_INDEPENDENT_CODE 1
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>)

#---Install shared build products---------------------------------------------------
install(TARGETS HepMC3protobufIO 
  DESTINATION ${CMAKE_INSTALL_LIBDIR} 
  COMPONENT protobufIOlibs)

install(FILES 
          ${CMAKE_CURRENT_SOURCE_DIR}/include/HepMC3/Writerprotobuf.h 
          ${CMAKE_CURRENT_SOURCE_DIR}/include/HepMC3/Readerprotobuf.h 
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/include/HepMC3 
  COMPONENT protobufIOdevel)

if (HEPMC3_BUILD_STATIC_LIBS)
#---Define static library---------------------------------------------------
  add_library(HepMC3protobufIO_static STATIC ${HepMC3protobufIO_sources})
  target_compile_definitions(HepMC3protobufIO_static PRIVATE HEPMC3_NO_EXPORTS)
  set_target_properties(HepMC3protobufIO_static 
    PROPERTIES 
      OUTPUT_NAME HepMC3protobufIO-static
      ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>
      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>)

#---Install static build products---------------------------------------------------
  install(TARGETS HepMC3protobufIO_static DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT protobufIOlibs)
endif()
