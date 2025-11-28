AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libhftrackefficiency.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --evelibs`

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

pkginclude_HEADERS = \
  HFTrackEfficiency.h

libhftrackefficiency_la_SOURCES = \
  HFTrackEfficiency.cc

libhftrackefficiency_la_LIBADD = \
  -lfun4all \
  -ldecayfinder_io \
  -ltrackbase_historic_io \
  -lphg4hit \
  -lphhepmc

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
	rm -f $(BUILT_SOURCES)
