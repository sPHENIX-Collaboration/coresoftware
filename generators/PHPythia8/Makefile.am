AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir` 

lib_LTLIBRARIES = libPHPythia8.la 

pkginclude_HEADERS = \
  PHPythia8.h \
  PHPy8GenTrigger.h \
  PHPy8ParticleTrigger.h \
  PHPy8JetTrigger.h 

libPHPythia8_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs` \
  `fastjet-config --libs`

libPHPythia8_la_LIBADD = \
  -lSubsysReco \
  -lpythia8 \
  -lphhepmc \
  -lHepMC

libPHPythia8_la_SOURCES = \
  PHPythia8.cc \
  PHPy8GenTrigger.cc \
  PHPy8ParticleTrigger.cc \
  PHPy8JetTrigger.cc

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHPythia8.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
