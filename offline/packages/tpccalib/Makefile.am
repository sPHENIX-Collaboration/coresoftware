##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
	libtpccalib.la

AM_CPPFLAGS = \
	-I$(includedir) \
	-I$(OFFLINE_MAIN)/include \
	-I$(ROOTSYS)/include \
	-I$(OFFLINE_MAIN)/include/eigen3

AM_LDFLAGS = \
	-L$(libdir) \
  -L$(ROOTSYS)/lib \
	-L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  	TpcSpaceChargeReconstruction.h

if ! MAKEROOT6
	ROOT5_DICTS = \
		TpcSpaceChargeReconstruction_Dict.cc
endif

libtpccalib_la_SOURCES = \
	$(ROOTDICTS) \
	TpcSpaceChargeReconstruction.cc

libtpccalib_la_LIBADD = \
	-lphool \
  -lSubsysReco \
  -lg4detectors_io \
	-ltrack_io \
  -ltrackbase_historic_io

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = \
	testexternals.cc

noinst_PROGRAMS = \
	testexternals_tpccalib

testexternals_tpccalib_SOURCES = testexternals.cc
testexternals_tpccalib_LDADD = libtpccalib.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
