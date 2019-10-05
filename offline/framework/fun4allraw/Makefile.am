AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

if ! MAKEROOT6
ROOT5_SUBSYSRECODICTS =

ROOT5_FUN4ALLDICTS = \
  Fun4AllPrdfInputManager_Dict.cc \
  Fun4AllPrdfOutputManager_Dict.cc \
  Fun4AllOutputManager_Dict.cc
endif

pkginclude_HEADERS = \
  Fun4AllEventOutStream.h \
  Fun4AllEventOutputManager.h \
  Fun4AllFileOutStream.h \
  Fun4AllPrdfInputManager.h \
  Fun4AllPrdfOutputManager.h \
  Fun4AllRolloverFileOutStream.h

lib_LTLIBRARIES = \
  libfun4all_raw.la

libfun4all_raw_la_SOURCES = \
  $(ROOT5_FUN4ALLDICTS) \
  Fun4AllEventOutStream.cc \
  Fun4AllEventOutputManager.cc \
  Fun4AllFileOutStream.cc \
  Fun4AllPrdfInputManager.cc \
  Fun4AllPrdfOutputManager.cc \
  Fun4AllRolloverFileOutStream.cc

libfun4all_raw_la_LIBADD = \
  -lfun4all \
  -lEvent \
  -lphool_raw

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libfun4all_raw.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
