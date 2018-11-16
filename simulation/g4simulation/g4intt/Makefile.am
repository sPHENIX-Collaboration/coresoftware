AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libg4intt_io.la \
    libg4intt.la

AM_CPPFLAGS = \
    -DCGAL_DISABLE_ROUNDING_MATH_CHECK \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(ROOTSYS)/include \
    -I$(G4_MAIN)/include \
    -I$(OPT_SPHENIX)/include

# set in configure.in to check if gcc version >= 4.8
#if GCC_GE_48
#  AM_CXXFLAGS = -std=c++11
#endif

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libg4intt_io_la_LIBADD = \
  -L$(ROOTSYS)/lib \
  -lCLHEP \
  -lphool \
  -lboost_filesystem \
  -lboost_system \
  -lpdbcalBase \
  -lXMLIO \
  -lphparameter \
  -lg4detectors 

libg4intt_la_LIBADD = \
  libg4intt_io.la \
  -lphg4gdml \
  -lphool  \
  -lCGAL \
  -lSubsysReco \
  -lg4testbench \
  -lg4detectors

pkginclude_HEADERS = \
   PHG4CylinderGeom_INTT.h \
  PHG4CylinderGeomINTT.h \
   PHG4INTTCellReco.h \
  PHG4INTTDefs.h \
  PHG4INTTSubsystem.h 
 
libg4intt_io_la_SOURCES = \
  PHG4CylinderGeom_INTT.cc \
  PHG4CylinderGeom_INTT_Dict.cc \
  PHG4CylinderGeomINTT.cc \
  PHG4CylinderGeomINTT_Dict.cc 

libg4intt_la_SOURCES = \
  PHG4INTTCellReco.cc \
  PHG4INTTCellReco_Dict.cc \
  PHG4INTTDefs_Dict.cc \
  PHG4INTTDetector.cc \
  PHG4INTTSteppingAction.cc \
  PHG4INTTSubsystem.cc \
  PHG4INTTSubsystem_Dict.cc \
  PHG4INTTParameterisation.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_g4intt \
  testexternals_g4intt_io

testexternals_g4intt_SOURCES = testexternals.cc
testexternals_g4intt_LDADD = libg4intt.la

testexternals_g4intt_io_SOURCES = testexternals.cc
testexternals_g4intt_io_LDADD = libg4intt_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
