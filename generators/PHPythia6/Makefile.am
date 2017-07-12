AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -I$(includedir) -I$(OFFLINE_MAIN)/include -I`root-config --incdir`

AM_CXXFLAGS = -Werror

lib_LTLIBRARIES = libPHPythia6.la

pkginclude_HEADERS = \
  PHPythia6.h \
  PHPy6GenTrigger.h \
  PHPy6ForwardElectronTrig.h \
  PHPy6ParticleTrigger.h

libPHPythia6_la_LDFLAGS =  \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  `root-config --libs`

libPHPythia6_la_LIBADD = \
  -lphool -lSubsysReco -lfun4all -lPythia6 -lEG -lEGPythia6 \
  -lphhepmc -lHepMC -lHepMCfio -lgsl -lgslcblas -lfastjet \
  -lCGAL


libPHPythia6_la_SOURCES = \
  PHPythia6.C \
  PHPythia6_Dict.C \
  PHPy6GenTrigger.C \
  PHPy6GenTrigger.h \
  PHPy6JetTrigger.C \
  PHPy6JetTrigger.h \
  PHPy6ForwardElectronTrig.C \
  PHPy6ForwardElectronTrig.h \
  PHPy6ParticleTrigger.C \
  PHPy6ParticleTrigger.h


BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHPythia6.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

PHPythia6_Dict.C: \
  PHPythia6.h \
  PHPy6GenTrigger.h \
  PHPy6JetTrigger.h \
  PHPy6ForwardElectronTrig.h \
  PHPy6ParticleTrigger.h \
  PHPythia6LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict*
	rm -f testexternals*

