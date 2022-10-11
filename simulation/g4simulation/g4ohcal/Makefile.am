AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

lib_LTLIBRARIES = \
  libg4ohcal.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libg4ohcal_la_LIBADD = \
  -lfun4all \
  -lphg4hit \
  -lg4detectors

pkginclude_HEADERS = \
  PHG4OHCalSubsystem.h

libg4ohcal_la_SOURCES = \
  PHG4OHCalDetector.cc \
  PHG4OHCalDisplayAction.cc \
  PHG4OHCalSteppingAction.cc \
  PHG4OHCalSubsystem.cc \
  PHG4OHCalFieldSetup.cc 

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libg4ohcal.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f $(BUILT_SOURCES)
