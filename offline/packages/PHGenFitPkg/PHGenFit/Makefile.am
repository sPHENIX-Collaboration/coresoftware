AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

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
  PlanarMeasurement.cc \
  SpacepointMeasurement.cc \
  Track.cc \
  Tools.cc

libPHGenFit_la_LIBADD = \
  -lgenfit2 \
  -lgenfit2exp \
  -lphfield

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHGenFit.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

