AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

lib_LTLIBRARIES = libvararray.la

pkginclude_HEADERS = \
  VariableArray.h \
  VariableArrayContainer.h \
  VariableArrayIds.h \
  VariableArrayUtils.h

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include 

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib

libvararray_la_LIBADD = \
  -lhalf \
  -lphool


libvararray_la_SOURCES = \
  VariableArray.cc \
  VariableArrayContainer.cc \
  VariableArrayDict.cc \
  VariableArrayUtils.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = \
  libvararray.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

VariableArrayDict.cc: \
  VariableArray.h \
  VariableArrayContainer.h \
  VariableArrayLinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^


clean-local:
	rm -f *Dict.* testexternals.cc *.pcm
