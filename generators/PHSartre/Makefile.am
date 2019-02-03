AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OPT_SPHENIX)/include \
  -I`root-config --incdir` 

AM_CXXFLAGS = -Werror

lib_LTLIBRARIES = libPHSartre.la 

pkginclude_HEADERS = \
  PHSartre.h \
  PHSartreGenTrigger.h \
  PHSartreParticleTrigger.h 

libPHSartre_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)lib \
  `root-config --libs`

libPHSartre_la_LIBADD = \
  -lsartre \
  -lgemini \
  -lSubsysReco \
  -lphhepmc \
  -lHepMC \
  -lUnuran \
  -lEG \
  -lMathMore

if ! MAKEROOT6
ROOT5_DICTS = \
  PHSartre_Dict.cc \
  PHSartreGenTrigger_Dict.cc \
  PHSartreParticleTrigger_Dict.cc
endif

libPHSartre_la_SOURCES = \
  $(ROOT5_DICTS) \
  PHSartre.cc \
  PHSartreGenTrigger.cc \
  PHSartreParticleTrigger.cc

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHSartre.la

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
