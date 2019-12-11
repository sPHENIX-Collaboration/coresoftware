AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = libPHPythia6.la

pkginclude_HEADERS = \
  PHPy6ForwardElectronTrig.h \
  PHPy6GenTrigger.h \
  PHPy6JetTrigger.h \
  PHPy6ParticleTrigger.h \
  PHPythia6.h

libPHPythia6_la_LDFLAGS =  \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  `root-config --libs`

libPHPythia6_la_LIBADD = \
  -lSubsysReco \
  -lEGPythia6 \
  -lphhepmc \
  -lHepMC \
  -lHepMCfio \
  -lfastjet \
  -lCGAL

if ! MAKEROOT6
ROOT5_DICTS = \
  PHPy6ForwardElectronTrig_Dict.cc \
  PHPy6GenTrigger_Dict.cc \
  PHPy6JetTrigger_Dict.cc \
  PHPy6ParticleTrigger_Dict.cc \
  PHPythia6_Dict.cc 
endif

libPHPythia6_la_SOURCES = \
  $(ROOT5_DICTS) \
  PHPythia6.cc \
  PHPy6GenTrigger.cc \
  PHPy6JetTrigger.cc \
  PHPy6ForwardElectronTrig.cc \
  PHPy6ParticleTrigger.cc


BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHPythia6.la

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

