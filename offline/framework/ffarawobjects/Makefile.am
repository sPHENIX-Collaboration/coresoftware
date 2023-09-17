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
  InttRawHit_Dict.cc \
  InttRawHitContainer_Dict.cc \
  InttRawHitContainerv1_Dict.cc \
  InttRawHitv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  InttRawHit_Dict_rdict.pcm \
  InttRawHitContainer_Dict_rdict.pcm \
  InttRawHitContainerv1_Dict_rdict.pcm \
  InttRawHitv1_Dict_rdict.pcm

pkginclude_HEADERS = \
  InttRawHit.h \
  InttRawHitContainer.h \
  InttRawHitContainerv1.h \
  InttRawHitv1.h

libffarawobjects_la_SOURCES = \
  $(ROOTDICTS) \
  InttRawHitContainerv1.cc

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
