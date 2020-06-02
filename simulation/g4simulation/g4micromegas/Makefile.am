AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include\
  -I$(G4_MAIN)/include 

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  PHG4MicromegasSubsystem.h

lib_LTLIBRARIES = \
  libPHG4Micromegas.la

libPHG4Micromegas_la_SOURCES = \
  PHG4MicromegasSubsystem.cc\
  PHG4MicromegasDetector.cc\
  PHG4MicromegasSteppingAction.cc

libPHG4Micromegas_la_LIBADD = \
  -lphool \
  -lSubsysReco\
  -lg4detectors\
  -lg4testbench 

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libPHG4Micromegas.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
