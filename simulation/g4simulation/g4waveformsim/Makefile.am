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
  CaloWaveformSim.h

lib_LTLIBRARIES = \
  libCaloWaveformSim.la

libCaloWaveformSim_la_SOURCES = \
  CaloWaveformSim.cc

libCaloWaveformSim_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lcalo_io \
  -lfun4all \
  -lg4detectors_io \
  -lg4detectors \
  -lsph_onnx \
  -lcalo_reco \
  -lcalo_io \
  -lcdbobjects \
  -lphg4hit

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libCaloWaveformSim.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
