AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I$(G4_MAIN)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  PHG4EPDetector.h \
  PHG4EPSteppingAction.h \
  PHG4EPDSubsystem.h

lib_LTLIBRARIES = \
  libg4epd.la

libg4epd_la_SOURCES = \
  PHG4EPDetector.cc \
  PHG4EPSteppingAction.cc \
  PHG4EPDSubsystem.cc

libg4epd_la_LIBADD = \
  -lphool \
  -lg4detectors

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternalsg4epd

testexternalsg4epd_SOURCES = testexternals.cc
testexternalsg4epd_LDADD = libg4epd.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
