AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libg4evtgendecayer.la

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
  -I`root-config --incdir` 


libg4evtgendecayer_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  `geant4-config --libs` \
  ../TestEvtGen/PHEvtGenDecayer.h

libg4evtgendecayer_la_LIBADD = \
  -lPythia6 \
  -lEvtGenExternal \
  -levtgendecayer

libg4evtgendecayer_la_SOURCES = \
  G4EvtGenDecayer.cc \
  EvtGenExtDecayerPhysics.cc

pkginclude_HEADERS = \
  EvtGenExtDecayerPhysics.hh

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libg4evtgendecayer.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
