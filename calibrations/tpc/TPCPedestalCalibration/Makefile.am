AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  TPCPedestalCalibration.h

lib_LTLIBRARIES = \
  libTPCPedestalCalibration.la

libTPCPedestalCalibration_la_SOURCES = \
  TPCPedestalCalibration.cc

libTPCPedestalCalibration_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lcdbobjects \
  -lsphenixnpc

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libTPCPedestalCalibration.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
