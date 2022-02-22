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
  PHG4MicromegasDigitizer.h \
  PHG4MicromegasHitReco.h \
  PHG4MicromegasSubsystem.h

lib_LTLIBRARIES = \
  libg4micromegas.la

libg4micromegas_la_SOURCES = \
  PHG4MicromegasDetector.cc \
  PHG4MicromegasDigitizer.cc \
  PHG4MicromegasDisplayAction.cc \
  PHG4MicromegasHitReco.cc \
  PHG4MicromegasSteppingAction.cc \
  PHG4MicromegasSubsystem.cc

libg4micromegas_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lg4detectors \
  -lg4testbench \
  -lmicromegas_io

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libg4micromegas.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
