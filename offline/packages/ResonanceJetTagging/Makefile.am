AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libresonancejettagging.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
	-DHomogeneousField

pkginclude_HEADERS = \
  ResonanceJetTagging.h

libresonancejettagging_la_SOURCES = \
  ResonanceJetTagging.cc

libresonancejettagging_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lcalo_io \
  -lfun4all \
  -lg4detectors_io \
  -lg4jets \
  -lkfparticle_sphenix_io \
  -lparticleflow_io \
  -lphg4hit \
  `fastjet-config --libs`


################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libresonancejettagging.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
