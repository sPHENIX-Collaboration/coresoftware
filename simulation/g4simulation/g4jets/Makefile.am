AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libg4jets.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4jets_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `fastjet-config --libs`

libg4jets_la_LIBADD = \
  -lfun4all \
  -ljetbase \
  -lphg4hit \
  -lcalo_io \
  -lRecursiveTools \
  -lphhepmc_io \
  -ltrackbase_historic_io

pkginclude_HEADERS = \
  TruthJetInput.h \
  JetHepMCLoader.h


libg4jets_la_SOURCES = \
  JetHepMCLoader.cc \
  TruthJetInput.cc

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4jets

BUILT_SOURCES = testexternals.cc

testexternals_g4jets_SOURCES = testexternals.cc
testexternals_g4jets_LDADD = libg4jets.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
