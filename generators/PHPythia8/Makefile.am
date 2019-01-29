AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir` 

AM_CXXFLAGS = -Werror

lib_LTLIBRARIES = libPHPythia8.la 

pkginclude_HEADERS = \
  PHPythia8.h \
  PHPy8GenTrigger.h \
  PHPy8ParticleTrigger.h \
  PHPy8JetTrigger.h 

libPHPythia8_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

libPHPythia8_la_LIBADD = \
  -lSubsysReco \
  -lpythia8 \
  -lphhepmc \
  -lHepMC \
  -lfastjet \
  -lCGAL

if ! MAKEROOT6
ROOT5_DICTS = \
  PHPythia8_Dict.cc \
  PHPy8GenTrigger_Dict.cc \
  PHPy8JetTrigger_Dict.cc \
  PHPy8ParticleTrigger_Dict.cc
endif

libPHPythia8_la_SOURCES = \
  $(ROOT5_DICTS) \
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

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^


clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
