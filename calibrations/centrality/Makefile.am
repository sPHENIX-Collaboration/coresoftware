AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

pkginclude_HEADERS = \
  CentralityCalibrator.h \
  MbdCentralityTree.h

lib_LTLIBRARIES = \
  libcentralitycalib.la

libcentralitycalib_la_SOURCES = \
  CentralityCalibrator.cc \
  MbdCentralityTree.cc

libcentralitycalib_la_LIBADD = \
  -lcdbobjects \
  -lffarawobjects \
  -lfun4all \
  -lglobalvertex_io \
  -lmbd_io \
  -lphool \
  -lSubsysReco \
  -lzdcinfo_io

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libcentralitycalib.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
