AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(srcdir) \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

pkginclude_HEADERS = \
  TpcConditions.h \
  TpcConditionsReco.h

lib_LTLIBRARIES = \
  libTpcConditions.la

libTpcConditions_la_LIBADD = \
  -lffamodules \
  -lffarawobjects \
  -lcdbobjects \
  -lphool \
  -lSubsysReco

libTpcConditions_la_SOURCES = \
  TpcConditionsReco.cc

################################################

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libTpcConditions.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
