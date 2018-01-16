AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libphparameter_io.la \
  libphparameter_util.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

libphparameter_io_la_SOURCES = \
  PHG4Parameters.cc \
  PHG4Parameters_Dict.cc \
  PHG4ParametersContainer.cc \
  PHG4ParametersContainer_Dict.cc 

libphparameter_util_la_SOURCES = \
  PHG4ParameterContainerInterface.cc \
  PHG4ParameterContainerInterface_Dict.cc \
  PHG4ParameterInterface.cc \
  PHG4ParameterInterface_Dict.cc 

libphparameter_io_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libphparameter_io_la_LIBADD = \
  -L$(ROOTSYS)/lib \
  -lphool \
  -lpdbcalBase \
  -lboost_filesystem \
  -lboost_system 

libphparameter_util_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libphparameter_util_la_LIBADD = \
  libphparameter_io.la 

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  PHG4ParameterContainerInterface.h \
  PHG4ParameterInterface.h \
  PHG4Parameters.h \
  PHG4ParametersContainer.h 
  
################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_phparameter_util \
  testexternals_phparameter_io


testexternals_phparameter_util_SOURCES = testexternals.C
testexternals_phparameter_util_LDADD = libphparameter_util.la

testexternals_phparameter_io_SOURCES = testexternals.C
testexternals_phparameter_io_LDADD = libphparameter_io.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
