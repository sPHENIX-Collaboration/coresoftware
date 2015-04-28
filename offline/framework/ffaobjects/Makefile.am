AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

lib_LTLIBRARIES = libffaobjects.la

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

AM_CXXFLAGS = -Wall -Werror
AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib

libffaobjects_la_LIBADD = \
  -lphool

pkginclude_HEADERS = \
  FlagSave.h \
  FlagSavev1.h \
  RunHeader.h \
  SyncDefs.h \
  SyncObject.h \
  SyncObjectv1.h \
  SyncObjectv2.h

noinst_HEADERS = \
  libffaobjects_LinkDef.h

libffaobjects_la_SOURCES = \
  FlagSave.cc \
  FlagSavev1.cc \
  RunHeader.cc \
  SyncObject.cc \
  SyncObjectv1.cc \
  SyncObjectv2.cc

nodist_libffaobjects_la_SOURCES = libffaobjects_Dict.cc

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

libffaobjects_Dict.cc: \
  FlagSave.h \
  FlagSavev1.h \
  RunHeader.h \
  SyncObject.h \
  SyncObjectv1.h \
  SyncObjectv2.h \
  libffaobjects_LinkDef.h
	rootcint -f $@ -c $(INCLUDES) $^

clean-local:
	rm -f *Dict.*
