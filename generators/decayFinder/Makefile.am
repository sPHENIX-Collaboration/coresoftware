AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libdecayfinder_io.la \
    libdecayfinder.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --evelibs`

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ROOTSYS)/include \
  -I$(OFFLINE_MAIN)/include

RINCLUDES = \
  -I$(includedir) \
  -I$(ROOTSYS)/include \
  -I$(OFFLINE_MAIN)/include

pkginclude_HEADERS = \
  DecayFinderContainer.h \
  DecayFinder.h

ROOTDICTS = \
  DecayFinderContainer_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  DecayFinderContainer_Dict_rdict.pcm

libdecayfinder_io_la_SOURCES = \
  $(ROOTDICTS) \
  DecayFinderContainer.cc

libdecayfinder_la_SOURCES = \
  DecayFinderContainer.cc \ 
  DecayFinder.cc

libdecayfinder_io_la_LIBADD = \
  -lphool

libdecayfinder_la_LDFLAGS = \
  libdecayfinder_io.la \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lfun4all \
  -lphhepmc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals \
  testexternals_io

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libdecayfinder.la

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD = libdecayfinder_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
