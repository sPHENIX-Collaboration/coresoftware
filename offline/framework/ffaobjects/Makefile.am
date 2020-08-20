AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libffaobjects.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libffaobjects_la_LIBADD = \
  -lphool

ROOTDICTS = \
  FlagSave_Dict.cc \
  FlagSavev1_Dict.cc \
  RunHeader_Dict.cc \
  SyncObject_Dict.cc \
  SyncObjectv1_Dict.cc \
  EventHeader_Dict.cc \
  EventHeaderv1_Dict.cc \
  EventHeaderv2_Dict.cc

# for root6 we need pcm and dictionaries but only for
# i/o classes. For root5 we need only dictionaries but
# those for i/o and classes available on the cmd line
# MAKEROOT6 is set in the configure.ac, true for root6
if MAKEROOT6
# this is a tweak to install files in $(libdir), automake refuses
# to install other files in libdir, this construct gets around this
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
endif

pkginclude_HEADERS = \
  FlagSave.h \
  FlagSavev1.h \
  RunHeader.h \
  SyncDefs.h \
  SyncObject.h \
  SyncObjectv1.h \
  EventHeader.h \
  EventHeaderv1.h \
  EventHeaderv2.h

libffaobjects_la_SOURCES = \
  $(ROOTDICTS) \
  FlagSavev1.cc \
  RunHeader.cc \
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
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f *Dict.* *.pcm $(BUILT_SOURCES)
