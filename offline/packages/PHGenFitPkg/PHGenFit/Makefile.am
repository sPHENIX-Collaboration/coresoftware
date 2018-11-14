AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir` \
  -I$(GENFIT)/include \
  -I$(OFFLINE_MAIN)/include/eigen3

libPHGenFit_la_LDFLAGS = \
  -nodefaultlibs \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs` \
  -L$(GENFIT)/lib

lib_LTLIBRARIES = \
  libPHGenFit.la
  
pkginclude_HEADERS = \
  Fitter.h \
  Measurement.h \
  PlanarMeasurement.h \
  SpacepointMeasurement.h \
  Track.h \
  Tools.h

libPHGenFit_la_SOURCES = \
  Fitter.cc \
  Measurement.cc \
  PlanarMeasurement.cc \
  SpacepointMeasurement.cc \
  Track.cc \
  Tools.cc

libPHGenFit_la_LIBADD = \
  -lgenfit2 \
  -lgenfit2exp \
  -lphfield

BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHGenFit.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

testexternals_SOURCES = testexternals.C
