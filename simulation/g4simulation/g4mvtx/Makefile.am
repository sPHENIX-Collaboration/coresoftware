AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

lib_LTLIBRARIES = \
  libg4mvtx.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libg4mvtx_la_LIBADD = \
  -lfun4all \
  -lphg4hit \
  -lg4detectors \
  -ltrack_io \
  -lmvtx_io

pkginclude_HEADERS = \
  PHG4MvtxDefs.h \
  PHG4MvtxHitReco.h \
  PHG4MvtxSubsystem.h \
  PHG4EICMvtxSubsystem.h \
  PHG4MvtxDigitizer.h

libg4mvtx_la_SOURCES = \
  PHG4MvtxHitReco.cc \
  PHG4MvtxDetector.cc \
  PHG4EICMvtxDetector.cc \
  PHG4MvtxDisplayAction.cc \
  PHG4MvtxSteppingAction.cc \
  PHG4EICMvtxSteppingAction.cc \
  PHG4MvtxSubsystem.cc \
  PHG4EICMvtxSubsystem.cc \
  PHG4MvtxDigitizer.cc \
  PHG4MvtxCable.cc \
  PHG4MvtxServiceStructure.cc \
  PHG4MvtxSupport.cc

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libg4mvtx.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f $(BUILT_SOURCES)
