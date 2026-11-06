AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(G4_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  PHG4ChimneySubsystem.h \
  PHG4ChimneyDetector.h \
  PHG4ChimneyDisplayAction.h \
  ChimneySetup.h

lib_LTLIBRARIES = \
  libg4chimney.la

libg4chimney_la_SOURCES = \
  PHG4ChimneySubsystem.cc \
  PHG4ChimneyDetector.cc \
  PHG4ChimneyDisplayAction.cc

libg4chimney_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lg4detectors \
  -lphparameter

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = testexternals
testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libg4chimney.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
