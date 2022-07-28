AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  fillSpaceChargeMaps.h \
  Shifter.h

lib_LTLIBRARIES = \
  libfillSpaceChargeMaps.la

libfillSpaceChargeMaps_la_SOURCES = \
  fillSpaceChargeMaps.cc \
  Shifter.cc

libfillSpaceChargeMaps_la_LIBADD = \
  -lfun4all \
  -lphool \
  -lphg4hit \
  -lg4detectors \
  -lg4testbench \
  -lSubsysReco

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libfillSpaceChargeMaps.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
