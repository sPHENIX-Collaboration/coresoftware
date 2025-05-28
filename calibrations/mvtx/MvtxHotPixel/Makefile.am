AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  MvtxFakeHitRate.h \
  MvtxHitMapWriter.h 
  
lib_LTLIBRARIES = \
  libmvtxcalib.la

libmvtxcalib_la_SOURCES = \
  MvtxFakeHitRate.cc \
  MvtxHitMapWriter.cc
  
libmvtxcalib_la_LIBADD = \
  -lphool \
  -ltrack_io \
  -lg4detectors \
  -ltrackbase_historic_io \
  -lmvtx \
  -lffarawobjects \
  -lcdbobjects \
  -lCLHEP \
  -lSubsysReco

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libmvtxcalib.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)

