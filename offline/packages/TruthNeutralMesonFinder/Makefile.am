AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libtruthneutralmeson_io.la \
    libtruthneutralmeson.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --evelibs`

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

pkginclude_HEADERS = \
  TruthNeutralMeson.h \
  TruthNeutralMesonv1.h \
  TruthNeutralMesonContainer.h \
  TruthNeutralMesonBuilder.h 

ROOTDICTS = \
  TruthNeutralMeson_Dict.cc \
  TruthNeutralMesonv1_Dict.cc \
  TruthNeutralMesonContainer_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  TruthNeutralMeson_Dict_rdict.pcm \
  TruthNeutralMesonv1_Dict_rdict.pcm \
  TruthNeutralMesonContainer_Dict_rdict.pcm

libtruthneutralmeson_io_la_SOURCES = \
  $(ROOTDICTS) \
  TruthNeutralMeson.cc \
  TruthNeutralMesonv1.cc \
  TruthNeutralMesonContainer.cc

libtruthneutralmeson_la_SOURCES = \
  TruthNeutralMesonBuilder.cc 

libtruthneutralmeson_io_la_LIBADD = \
  -lphool

libtruthneutralmeson_la_LIBADD = \
  libtruthneutralmeson_io.la \
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
testexternals_LDADD = libtruthneutralmeson.la

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD = libtruthneutralmeson_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
