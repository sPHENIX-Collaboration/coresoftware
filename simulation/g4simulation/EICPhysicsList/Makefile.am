AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libEICPhysicsList.la

# set in configure.in to check if gcc version >= 4.8
# leave this here to show how it can be done, -std=c++11 is now
# enabled by default in our config.site
#if GCC_GE_48
#  AM_CXXFLAGS = -std=c++11
#endif

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I${G4_MAIN}/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs`

#libEICPhysicsList_la_LIBADD = \
#  -lPythia6

libEICPhysicsList_la_SOURCES = \
  AntiBaryonPhysics.cc \
  EICPhysicsList.cc \
  GammaLeptoNuclearPhysics.cc \
  HyperonPhysics.cc \
  IonPhysics.cc \
  KaonPhysics.cc \
  NeutronPhysics.cc \
  PionPhysics.cc \
  ProtonPhysics.cc

pkginclude_HEADERS = \
  EICPhysicsList.hh

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libEICPhysicsList.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
