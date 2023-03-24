AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libg4decayer.la

# set in configure.in to check if gcc version >= 4.8
# leave this here to show how it can be done, -std=c++11 is now
# enabled by default in our config.site
#if GCC_GE_48
#  AM_CXXFLAGS = -std=c++11
#endif

AM_CXXFLAGS=`geant4-config --cflags`

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -DEVTGEN_HEPMC3

libg4decayer_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `geant4-config --libs`

libg4decayer_la_LIBADD = \
  -lPythia6 \
  -lEvtGen \
  -lEvtGenExternal \
  -lHepMC3 \
  -lgsl \
  -lgslcblas \
  -lfun4all 


libg4decayer_la_SOURCES = \
  pythia6_common_address.c \
  Pythia6.cc \
  G4Pythia6Decayer.cc \
  G4Pythia6DecayerMessenger.cc \
  P6DExtDecayerPhysics.cc \
  G4EvtGenDecayer.cc \
  EvtGenExtDecayerPhysics.cc \
  PHEvtGenRandomEngine.cc 

pkginclude_HEADERS = \
  EDecayType.hh \
  P6DExtDecayerPhysics.hh \
  G4EvtGenDecayer.hh \
  EvtGenExtDecayerPhysics.hh \
  PHEvtGenRandomEngine.hh 

################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libg4decayer.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
