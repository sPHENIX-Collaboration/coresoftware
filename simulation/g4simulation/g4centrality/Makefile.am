AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libg4centrality.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4centrality_la_LIBADD = \
  -lcentrality_io \
  -lfun4all \
  -lphg4hit \
  -lphparameter

pkginclude_HEADERS = \
  PHG4CentralityReco.h

libg4centrality_la_SOURCES = \
  PHG4CentralityReco.cc

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4centrality

BUILT_SOURCES = testexternals.cc

testexternals_g4centrality_SOURCES = testexternals.cc
testexternals_g4centrality_LDADD = libg4centrality.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
