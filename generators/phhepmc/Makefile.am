AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
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
  HepMCFlowAfterBurner.h \
  PHGenIntegral.h \
  PHGenIntegralv1.h \
  PHHepMCDefs.h \
  PHHepMCGenEvent.h \
  PHHepMCGenEventv1.h \
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
  -lflowafterburner \
  -lgsl \
  -lgslcblas

ROOTDICTS = \
  PHGenIntegral_Dict.cc \
  PHGenIntegralv1_Dict.cc \
  PHHepMC_io_Dict.cc \
  PHHepMCGenEvent_Dict.cc \
  PHHepMCGenEventv1_Dict.cc \
  PHHepMCGenEventMap_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libphhepmc_la_SOURCES = \
  Fun4AllHepMCInputManager.cc \
  Fun4AllHepMCPileupInputManager.cc \
  Fun4AllHepMCOutputManager.cc \
  Fun4AllOscarInputManager.cc \
  HepMCFlowAfterBurner.cc \
  PHHepMCGenHelper.cc \
  PHHepMCParticleSelectorDecayProductChain.cc

libphhepmc_io_la_SOURCES = \
  $(ROOTDICTS) \
  PHHepMCGenEvent.cc \
  PHHepMCGenEventv1.cc \
  PHHepMCGenEventMap.cc \
  PHGenIntegral.cc \
  PHGenIntegralv1.cc


libphhepmc_io_la_LDFLAGS = \
  ${AM_LDFLAGS} \
  `root-config --libs`

libphhepmc_io_la_LIBADD = \
  -lphool \
  -lHepMC \
  -lCLHEP

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

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

PHHepMC_io_Dict.cc: \
  $(OFFLINE_MAIN)/include/HepMC/*.h \
  PHHepMC_ioLinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm 
