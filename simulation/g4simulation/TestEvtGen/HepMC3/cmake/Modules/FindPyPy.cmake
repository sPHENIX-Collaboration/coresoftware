find_program(PyPy_EXE  NAMES pypy${PyPy_VERSION} )

execute_process( COMMAND ${PyPy_EXE} -c "import  sys;print(sys.prefix)"  OUTPUT_VARIABLE pypy_prefix OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process( COMMAND ${PyPy_EXE} -c "import  distutils;print(distutils.sysconfig.get_python_inc())"  OUTPUT_VARIABLE pypy_inc OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(  COMMAND ${PyPy_EXE} -c "import  sys;v=sys.version_info;print(str(v[0])+'.'+str(v[1])+'.'+str(v[2]))"  OUTPUT_VARIABLE pypy_python_version OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process( COMMAND ${PyPy_EXE} -c "import  sys;v=sys.pypy_version_info;print(str(v[0])+'.'+str(v[1])+'.'+str(v[2]))" OUTPUT_VARIABLE PyPy_PyPy_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

set( PyPy_Python_VERSION ${pypy_python_version})
string(REPLACE "." ";" MAJMINMICRO "${PyPy_Python_VERSION}")
list(GET MAJMINMICRO 0 PyPy_Python_VERSION_MAJOR)
list(GET MAJMINMICRO 1 PyPy_Python_VERSION_MINOR)
list(GET MAJMINMICRO 2 PyPy_Python_VERSION_MICRO)

string(REPLACE "." ";" MAJMINMICRO2 "${PyPy_PyPy_VERSION}")
list(GET MAJMINMICRO2 0 PyPy_PyPy_VERSION_MAJOR)
list(GET MAJMINMICRO2 1 PyPy_PyPy_VERSION_MINOR)
list(GET MAJMINMICRO2 2 PyPy_PyPy_VERSION_MICRO)


find_library(PyPy_LIBRARY NAMES pypy-c )
set( PyPy_INCLUDE_DIR  ${pypy_inc})
set( PyPy_SITEARCH  ${pypy_prefix}/site-packages)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PyPy REQUIRED_VARS PyPy_INCLUDE_DIR PyPy_LIBRARY)
set(PyPy_INCLUDE_DIRS ${PyPy_INCLUDE_DIR})
set(PyPy_LIBRARIES ${PyPy_LIBRARY})
mark_as_advanced(PyPy_INCLUDE_DIR PyPy_INCLUDE_DIRS PyPy_LIBRARY PyPy_LIBRARIES PyPy_Python_VERSION_MAJOR)
