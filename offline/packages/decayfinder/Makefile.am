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
  -isystem$(ROOTSYS)/include \
  -I$(OFFLINE_MAIN)/include

pkginclude_HEADERS = \
  DecayFinderContainerBase.h \
  DecayFinderContainer_v1.h \
  DecayFinder.h

ROOTDICTS = \
  DecayFinderContainerBase_Dict.cc \
  DecayFinderContainer_v1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  DecayFinderContainerBase_Dict_rdict.pcm \
  DecayFinderContainer_v1_Dict_rdict.pcm

libdecayfinder_io_la_SOURCES = \
  $(ROOTDICTS) \
  DecayFinderContainerBase.cc \
  DecayFinderContainer_v1.cc 

libdecayfinder_la_SOURCES = \
  DecayFinder.cc

libdecayfinder_io_la_LIBADD = \
  -lphool

libdecayfinder_la_LIBADD = \
  libdecayfinder_io.la \
  -lfun4all \
  -lphhepmc \
  -lphg4hit

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
