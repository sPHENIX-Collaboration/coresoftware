AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  dEdxFitter.h \
  GlobaldEdxFitter.h \
  bethe_bloch.h

lib_LTLIBRARIES = \
  libdedxfitter.la

libdedxfitter_la_SOURCES = \
  dEdxFitter.cc \
  GlobaldEdxFitter.cc

libdedxfitter_la_LIBADD = \
  -lphool \
  -ltrack_io \
  -lg4detectors \
  -ltrackbase_historic_io \
  -ltrack_reco \
  -lglobalvertex \
  -lSubsysReco

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libdedxfitter.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
