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
  TpcSpaceChargeReconstruction.h \
  TpcSpaceChargeReconstructionHelper.h

libtpccalib_la_SOURCES = \
  TpcSpaceChargeReconstruction.cc \
  TpcSpaceChargeReconstructionHelper.cc

libtpccalib_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lg4detectors_io \
  -ltrack_io \
  -ltrackbase_historic_io

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
	rm -f $(BUILT_SOURCES)
