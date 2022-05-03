AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  PHG4EPDDetector.h \
  PHG4EPDSteppingAction.h \
  PHG4EPDSubsystem.h

lib_LTLIBRARIES = \
  libg4epd.la

libg4epd_la_SOURCES = \
  PHG4EPDDetector.cc \
  PHG4EPDDisplayAction.cc \
  PHG4EPDSteppingAction.cc \
  PHG4EPDSubsystem.cc

libg4epd_la_LIBADD = \
  -lphool \
  -lg4detectors

# leave generic name for test executable, they stay in the local build area
# this preserves a similar "look and feel"
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libg4epd.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
