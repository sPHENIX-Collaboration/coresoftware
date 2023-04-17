##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# list of shared libraries to produce
lib_LTLIBRARIES = \
  libepd_io.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(ROOTSYS)/lib

libepd_io_la_LIBADD = \
  -lphool \
  -lcalo_io

pkginclude_HEADERS = \
  EPDDefs.h \
  EpdGeom.h \
  EpdGeomV1.h

ROOTDICTS = \
  EpdGeom_Dict.cc \
  EpdGeomV1_Dict.cc

pcmdir = $(libdir)

nobase_dist_pcm_DATA = \
  EpdGeom_Dict_rdict.pcm \
  EpdGeomV1_Dict_rdict.pcm 

libepd_io_la_SOURCES = \
  $(ROOTDICTS) \
  EpdGeomV1.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_io

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD = libepd_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f $(BUILT_SOURCES)
