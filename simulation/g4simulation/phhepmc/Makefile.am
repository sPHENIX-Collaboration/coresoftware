AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`
#  `geant4-config --cflags-without-gui`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib

lib_LTLIBRARIES = libphhepmc.la

pkginclude_HEADERS = \
  Fun4AllHepMCInputManager.h \
  Fun4AllHepMCOutputManager.h \
  Fun4AllOscarInputManager.h \
  PHGenEvent.h \
  PHGenEventv1.h \
  PHGenEventList.h \
  PHGenEventListv1.h \
  PHHepMCGenEvent.h

libphhepmc_la_LDFLAGS = ${AM_LDFLAGS} `root-config --libs`
libphhepmc_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lHepMC \
  -lboost_iostreams \
  -lfun4all

libphhepmc_la_SOURCES = \
  Fun4AllHepMCInputManager.cc \
  Fun4AllHepMCOutputManager.cc \
  Fun4AllOscarInputManager.cc \
  PHHepMCParticleSelectorDecayProductChain.cc \
  PHGenEventv1.cc \
  PHGenEventListv1.cc \
  PHHepMCGenEvent.cc \
  PHHepMC_Dict.cc

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libphhepmc.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

PHHepMC_Dict.cc: \
  Fun4AllHepMCInputManager.h \
  Fun4AllHepMCOutputManager.h \
  Fun4AllOscarInputManager.h \
  HepMCFlowAfterBurner.h \
  PHHepMCParticleSelectorDecayProductChain.h \
  PHGenEvent.h \
  PHGenEventv1.h \
  PHGenEventList.h \
  PHGenEventListv1.h \
  PHHepMCGenEvent.h \
  $(OFFLINE_MAIN)/include/HepMC/*.h \
  PHHepMCLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict*
	rm -f testexternals*
