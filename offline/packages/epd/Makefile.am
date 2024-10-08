AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libepd_io.la \
   libepd.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libepd_io_la_LIBADD = \
  -lphool \
  -lcalo_io

libepd_la_LIBADD = \
  libepd_io.la \
  -lcdbobjects \
  -lfun4all \
  -lglobalvertex_io 
 
pkginclude_HEADERS = \
  EPDDefs.h \
  EpdGeom.h \
  EpdGeomV1.h \
  EpdReco.h

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


libepd_la_SOURCES = \
  EpdReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_epd_io \
  testexternals_epd

BUILT_SOURCES = \
  testexternals.cc

testexternals_epd_io_SOURCES = testexternals.cc
testexternals_epd_io_LDADD = libepd_io.la

testexternals_epd_SOURCES = testexternals.cc
testexternals_epd_LDADD = libepd.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
