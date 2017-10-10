AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libg4tpc.la

AM_CPPFLAGS = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(G4_MAIN)/include \
    -I$(ROOTSYS)/include

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libg4tpc_la_LIBADD = \
  -lphool \
  -lg4testbench \
  -lphg4hit \
  -lg4detectors

#pkginclude_HEADERS =

libg4tpc_la_SOURCES = \
  PHG4TPCDetector.cc \
  PHG4TPCSteppingAction.cc \
  PHG4TPCSubsystem.cc \
  PHG4TPCSubsystem_Dict.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^


################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libg4tpc.la

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
