AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

lib_LTLIBRARIES = \
  libeASTPhysicsList.la

# set in configure.in to check if gcc version >= 4.8
# leave this here to show how it can be done, -std=c++11 is now
# enabled by default in our config.site
#if GCC_GE_48
#  AM_CXXFLAGS = -std=c++11
#endif

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs`

libeASTPhysicsList_la_SOURCES = \
  eASTAntiBaryonPhysics.cc \
  eASTGammaLeptoNuclearPhysics.cc \
  eASTHyperonPhysics.cc \
  eASTIonPhysics.cc \
  eASTKaonPhysics.cc \
  eASTNeutronPhysics.cc \
  eASTPhysicsList.cc \
  eASTPhysicsListMessenger.cc \
  eASTPionPhysics.cc \
  eASTProtonPhysics.cc

pkginclude_HEADERS = \
  eASTPhysicsList.hh

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libeASTPhysicsList.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
