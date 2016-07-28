AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir` 

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib 

AM_CXXFLAGS = -Wall -Werror

lib_LTLIBRARIES = libPHPythia8.la 

pkginclude_HEADERS = \
  PHPythia8.h \
  PHPy8GenTrigger.h \
  PHPy8ParticleTrigger.h \
  PHPy8JetTrigger.h 

libPHPythia8_la_LDFLAGS = \
  ${AM_LDFLAGS} \
  `root-config --libs`

libPHPythia8_la_LIBADD = \
  -lgsl \
  -lgslcblas \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lboost_iostreams \
  -lpythia8 \
  -lphhepmc \
  -lHepMC \
  -lfastjet \
  -lCGAL

libPHPythia8_la_SOURCES = \
  PHPythia8.C \
  PHPy8GenTrigger.C \
  PHPy8ParticleTrigger.C \
  PHPy8JetTrigger.C \
  PHPythia8_Dict.C

BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHPythia8.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

PHPythia8_Dict.C: \
  PHPythia8.h \
  PHPy8GenTrigger.h \
  PHPy8ParticleTrigger.h \
  PHPy8JetTrigger.h \
  PHPythia8LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict*
	rm -f testexternals*
