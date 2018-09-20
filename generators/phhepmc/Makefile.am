AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I/opt/sphenix/core/HepMC-3.0.0/include \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L/opt/sphenix/core/HepMC-3.0.0/lib64 \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib

lib_LTLIBRARIES = \
  libphhepmc_io.la \
  libphhepmc.la

pkginclude_HEADERS = \
  Fun4AllHepMCInputManager.h \
  Fun4AllHepMCPileupInputManager.h \
  Fun4AllHepMCOutputManager.h \
  Fun4AllOscarInputManager.h \
  PHGenEvent.h \
  PHGenEventv1.h \
  PHGenEventList.h \
  PHGenEventListv1.h \
  PHGenIntegral.h \
  PHGenIntegralv1.h \
  PHHepMCGenEvent.h \
  PHHepMCGenEventMap.h \
  PHHepMCGenHelper.h

libphhepmc_la_LDFLAGS = \
  ${AM_LDFLAGS} \
  `root-config --libs`

libphhepmc_la_LIBADD = \
  libphhepmc_io.la \
  -lphool \
  -lSubsysReco \
  -lboost_iostreams \
  -lfun4all \
  -lgsl \
  -lgslcblas

libphhepmc_la_SOURCES = \
  Fun4AllHepMCInputManager.cc \
  Fun4AllHepMCInputManager_Dict.C \
  Fun4AllHepMCPileupInputManager.cc \
  Fun4AllHepMCPileupInputManager_Dict.C \
  Fun4AllHepMCOutputManager.cc \
  Fun4AllHepMCOutputManager_Dict.C \
  Fun4AllOscarInputManager.cc \
  Fun4AllOscarInputManager_Dict.C \
  PHHepMCGenHelper.cc \
  PHHepMCGenHelper_Dict.C

libphhepmc_io_la_SOURCES = \
  PHGenEventv1.cc \
  PHGenEventListv1.cc \
  PHHepMCGenEvent.cc \
  PHHepMCGenEvent_Dict.C \
  PHHepMCGenEventMap.cc \
  PHHepMCGenEventMap_Dict.C \
  PHGenIntegral.cc \
  PHGenIntegralv1.cc

libphhepmc_io_la_LDFLAGS = \
  ${AM_LDFLAGS} \
  `root-config --libs`

libphhepmc_io_la_LIBADD = \
  -lphool \
  -lHepMC

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_libphhepmc \
  testexternals_libphhepmc_io

testexternals_libphhepmc_SOURCES = testexternals.cc
testexternals_libphhepmc_LDADD = libphhepmc.la

testexternals_libphhepmc_io_SOURCES = testexternals.cc
testexternals_libphhepmc_io_LDADD = libphhepmc_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(AM_CPPFLAGS) $^

PHHepMC_io_Dict.C: \
  PHGenEvent.h \
  PHGenEventv1.h \
  PHGenEventList.h \
  PHGenEventListv1.h \
  PHHepMCGenEvent.h \
  PHHepMCGenEventMap.h \
  /opt/sphenix/core/HepMC-3.0.0/include/HepMC/*.h \
  PHHepMCLinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(AM_CPPFLAGS) $^



clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
