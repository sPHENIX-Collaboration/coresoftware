AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libvararray.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include 

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libvararray_la_LIBADD = \
  -lhalf \
  -lphool

pkginclude_HEADERS = \
  VariableArray.h \
  VariableArrayContainer.h \
  VariableArrayIds.h \
  VariableArrayUtils.h

ROOTDICTS = \
  VariableArray_Dict.cc \
  VariableArrayContainer_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  VariableArray_Dict_rdict.pcm \
  VariableArrayContainer_Dict_rdict.pcm

libvararray_la_SOURCES = \
  $(ROOTDICTS) \
  VariableArray.cc \
  VariableArrayContainer.cc \
  VariableArrayUtils.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking test to make sure we do not have unresolved symbols

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libvararray.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict.* testexternals.cc *.pcm
