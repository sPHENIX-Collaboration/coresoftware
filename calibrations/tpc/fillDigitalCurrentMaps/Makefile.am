AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  readDigitalCurrents.h

lib_LTLIBRARIES = \
  libreadDigitalCurrents.la

libreadDigitalCurrents_la_SOURCES = \
  readDigitalCurrents.cc

libreadDigitalCurrents_la_LIBADD = \
  -lphool \
  -lfun4all \
  -lphg4hit \
  -ltrackbase_historic_io \
  -ltrack_io \
  -ltpc_io \
  -lg4detectors \
  -lg4testbench \
  -lSubsysReco

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libreadDigitalCurrents.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
