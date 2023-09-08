AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libffarawobjects.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libffarawobjects_la_LIBADD = \
  -lphool

ROOTDICTS = \
  InttHit_Dict.cc \
  InttHitContainer_Dict.cc \
  InttHitContainerv1_Dict.cc \
  InttHitv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  InttHit_Dict_rdict.pcm \
  InttHitContainer_Dict_rdict.pcm \
  InttHitContainerv1_Dict_rdict.pcm \
  InttHitv1_Dict_rdict.pcm

pkginclude_HEADERS = \
  InttHit.h \
  InttHitContainer.h \
  InttHitContainerv1.h \
  InttHitv1.h

libffarawobjects_la_SOURCES = \
  $(ROOTDICTS) \
  InttHitContainerv1.cc

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libffarawobjects.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f *Dict.* *.pcm $(BUILT_SOURCES)
