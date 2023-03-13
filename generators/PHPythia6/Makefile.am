AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

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
  `root-config --libs` \
  `fastjet-config --libs`

libPHPythia6_la_LIBADD = \
  -lSubsysReco \
  -lEGPythia6 \
  -lphhepmc \
  -lHepMC \
  -lHepMCfio

libPHPythia6_la_SOURCES = \
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

clean-local:
	rm -f $(BUILT_SOURCES)

