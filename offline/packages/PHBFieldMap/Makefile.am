AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.C

lib_LTLIBRARIES = libPHBFieldMap.la

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I`root-config --incdir`

libPHBFieldMap_la_LDFLAGS = \
  `root-config --glibs`

noinst_HEADERS = \
  PHBFieldMapLinkDef.h

include_HEADERS = \
  PHBFieldMap.hh \
  PHBFieldMapObj.h

libPHBFieldMap_la_SOURCES = \
  PHBFieldMap.C \
  PHBFieldMapObj.cc \
  PHBFieldMapDict.C

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libPHBFieldMap.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

%Dict.C: %.hh %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
