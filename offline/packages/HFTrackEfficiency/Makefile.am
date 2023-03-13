AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libhftrackefficiency.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --evelibs`

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(ROOTSYS)/include \
  -I$(OFFLINE_MAIN)/include

pkginclude_HEADERS = \
  HFTrackEfficiency.h

pcmdir = $(libdir)

libhftrackefficiency_la_SOURCES = \
  HFTrackEfficiency.cc

libhftrackefficiency_la_LIBADD = \
  -lfun4all \
  -ldecayfinder_io \
  -ltrackbase_historic_io \
  -lphg4hit \
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
testexternals_LDADD = libhftrackefficiency.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
