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
  -lgsl \
  -lgslcblas \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lboost_iostreams \
  -lphhepmc \
  -lHepMC \
  -lUnuran \
  -lEG \
  -lMathMore

libPHSartre_la_SOURCES = \
  PHSartre.C \
  PHSartreGenTrigger.C \
  PHSartreParticleTrigger.C \
  PHSartre_Dict.C

BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHSartre.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

PHSartre_Dict.C: \
  PHSartre.h \
  PHSartreGenTrigger.h \
  PHSartreParticleTrigger.h \
  PHSartreLinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict*
	rm -f testexternals*
