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
  InttRawHitv1_Dict.cc \
  MicromegasRawHit_Dict.cc \
  MicromegasRawHitContainer_Dict.cc \
  MicromegasRawHitContainerv1_Dict.cc \
  MicromegasRawHitv1_Dict.cc \
  MvtxRawHit_Dict.cc \
  MvtxRawHitContainer_Dict.cc \
  MvtxRawHitContainerv1_Dict.cc \
  MvtxRawHitv1_Dict.cc \
  MvtxRawEvtHeader_Dict.cc \
  MvtxRawEvtHeaderv1_Dict.cc \
  TpcRawHit_Dict.cc \
  TpcRawHitContainer_Dict.cc \
  TpcRawHitContainerv1_Dict.cc \
  TpcRawHitv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  InttRawHit_Dict_rdict.pcm \
  InttRawHitContainer_Dict_rdict.pcm \
  InttRawHitContainerv1_Dict_rdict.pcm \
  InttRawHitv1_Dict_rdict.pcm \
  MicromegasRawHit_Dict_rdict.pcm \
  MicromegasRawHitContainer_Dict_rdict.pcm \
  MicromegasRawHitContainerv1_Dict_rdict.pcm \
  MicromegasRawHitv1_Dict_rdict.pcm \
  MvtxRawHit_Dict_rdict.pcm \
  MvtxRawHitv1_Dict_rdict.pcm \
  MvtxRawHitContainer_Dict_rdict.pcm \
  MvtxRawHitContainerv1_Dict_rdict.pcm \
  MvtxRawEvtHeader_Dict_rdict.pcm \
  MvtxRawEvtHeaderv1_Dict_rdict.pcm \
  TpcRawHit_Dict_rdict.pcm \
  TpcRawHitContainer_Dict_rdict.pcm \
  TpcRawHitContainerv1_Dict_rdict.pcm \
  TpcRawHitv1_Dict_rdict.pcm

pkginclude_HEADERS = \
  InttRawHit.h \
  InttRawHitContainer.h \
  InttRawHitContainerv1.h \
  InttRawHitv1.h \
  MicromegasRawHit.h \
  MicromegasRawHitContainer.h \
  MicromegasRawHitContainerv1.h \
  MicromegasRawHitv1.h \
  MvtxRawHit.h \
  MvtxRawHitv1.h \
  MvtxRawHitContainer.h \
  MvtxRawHitContainerv1.h \
  MvtxRawEvtHeader.h \
  MvtxRawEvtHeaderv1.h \
  TpcRawHit.h \
  TpcRawHitContainer.h \
  TpcRawHitContainerv1.h \
  TpcRawHitv1.h

libffarawobjects_la_SOURCES = \
  $(ROOTDICTS) \
  InttRawHitContainerv1.cc \
  InttRawHitv1.cc \
  MicromegasRawHitContainerv1.cc \
  MicromegasRawHitv1.cc \
  MvtxRawHitContainerv1.cc \
  MvtxRawHitv1.cc \
  MvtxRawEvtHeaderv1.cc \
  TpcRawHitContainerv1.cc \
  TpcRawHitv1.cc

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
