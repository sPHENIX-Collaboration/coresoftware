AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libdecayfinder.la

AM_LDFLAGS = \
  -L$(libdir) \
  `root-config --libs`

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ROOTSYS)/include \
  -I$(OFFLINE_MAIN)/include 

pkginclude_HEADERS = \
  DecayFinder.h

libdecayfinder_la_SOURCES = \
  DecayFinder.cc

libdecayfinder_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -I$(ROOTSYS)/include \
  -lfun4all \
  -lg4eval \
  -lkfparticle_sphenix \
  -lphhepmc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals 

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libdecayfinder.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
