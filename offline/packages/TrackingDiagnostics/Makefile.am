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
  KshortReconstruction.h \
  helixResiduals.h \
  TrackResiduals.h \
  TrackContainerCombiner.h \
  TrackSeedTrackMapConverter.h \
  TrkrNtuplizer.h

lib_LTLIBRARIES = \
  libTrackingDiagnostics.la

libTrackingDiagnostics_la_SOURCES = \
  KshortReconstruction.cc \
  helixResiduals.cc \
  TrackResiduals.cc \
  TrackContainerCombiner.cc \
  TrackSeedTrackMapConverter.cc \
  TrkrNtuplizer.cc

libTrackingDiagnostics_la_LIBADD = \
  -lphool \
  -lfun4all \
  -lg4eval \
  -lg4eval_io \
  -lphg4hit \
  -lg4detectors_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -ltpc_io \
  -lmvtx_io \
  -lintt_io \
  -lmicromegas_io \
  -lSubsysReco \
  -ltrackeralign

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libTrackingDiagnostics.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
