AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libffaobjects_io.la \
  libffaobjects.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libffaobjects_io_la_LIBADD = \
  -lphool

libffaobjects_la_LIBADD = \
  libffaobjects_io.la \
  -lfun4all \
  -lSubsysReco

ROOTDICTS = \
  FlagSave_Dict.cc \
  FlagSavev1_Dict.cc \
  RunHeader_Dict.cc \
  SyncObject_Dict.cc \
  SyncObjectv1_Dict.cc \
  EventHeader_Dict.cc \
  EventHeaderv1_Dict.cc \
  EventHeaderv2_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  FlagSave_Dict_rdict.pcm \
  FlagSavev1_Dict_rdict.pcm \
  RunHeader_Dict_rdict.pcm \
  SyncObject_Dict_rdict.pcm \
  SyncObjectv1_Dict_rdict.pcm \
  EventHeader_Dict_rdict.pcm \
  EventHeaderv1_Dict_rdict.pcm \
  EventHeaderv2_Dict_rdict.pcm

pkginclude_HEADERS = \
  FlagSave.h \
  FlagSavev1.h \
  RunHeader.h \
  SyncDefs.h \
  SyncObject.h \
  SyncObjectv1.h \
  EventHeader.h \
  EventHeaderv1.h \
  EventHeaderv2.h \
  SyncReco.h

libffaobjects_io_la_SOURCES = \
  $(ROOTDICTS) \
  FlagSavev1.cc \
  RunHeader.cc \
  SyncObject.cc \
  SyncObjectv1.cc \
  EventHeader.cc \
  EventHeaderv1.cc \
  EventHeaderv2.cc

libffaobjects_la_SOURCES = \
  SyncReco.cc

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_io \
  testexternals

testexternals_io_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

testexternals_io_LDADD = \
  libffaobjects_io.la

testexternals_LDADD = \
  libffaobjects.la

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
