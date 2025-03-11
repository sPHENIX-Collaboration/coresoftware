AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libffaobjects.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libffaobjects_la_LIBADD = \
  -lphool

ROOTDICTS = \
  CdbUrlSave_Dict.cc \
  CdbUrlSavev1_Dict.cc \
  FlagSave_Dict.cc \
  FlagSavev1_Dict.cc \
  RunHeader_Dict.cc \
  RunHeaderv1_Dict.cc \
  SyncObject_Dict.cc \
  SyncObjectv1_Dict.cc \
  EventHeader_Dict.cc \
  EventHeaderv1_Dict.cc \
  EventHeaderv2_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

pkginclude_HEADERS = \
  CdbUrlSave.h \
  CdbUrlSavev1.h \
  FlagSave.h \
  FlagSavev1.h \
  RunHeader.h \
  RunHeaderv1.h \
  SyncDefs.h \
  SyncObject.h \
  SyncObjectv1.h \
  EventHeader.h \
  EventHeaderv1.h \
  EventHeaderv2.h

libffaobjects_la_SOURCES = \
  $(ROOTDICTS) \
  CdbUrlSave.cc \
  CdbUrlSavev1.cc \
  FlagSavev1.cc \
  RunHeader.cc \
  RunHeaderv1.cc \
  SyncObject.cc \
  SyncObjectv1.cc \
  EventHeader.cc \
  EventHeaderv1.cc \
  EventHeaderv2.cc

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

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
