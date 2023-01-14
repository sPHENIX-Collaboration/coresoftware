AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
   libphgeom_io.la \
   libphgeom.la

pkginclude_HEADERS = \
  PHGeomFileImport.h \
  PHGeomIOTGeo.h \
  PHGeomTGeo.h \
  PHGeomUtility.h

ROOTDICTS = \
  PHGeomIOTGeo_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
    PHGeomIOTGeo_Dict_rdict.pcm

libphgeom_la_SOURCES = \
  PHGeomTGeo.cc \
  PHGeomFileImport.cc \
  PHGeomUtility.cc

libphgeom_io_la_SOURCES = \
  $(ROOTDICTS) \
  PHGeomIOTGeo.cc

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libphgeom_la_LIBADD = \
  libphgeom_io.la \
  -lSubsysReco \
  -luuid

libphgeom_io_la_LIBADD = \
  -lphool \
  -lGeom

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_libphgeom \
  testexternals_libphgeom_io

testexternals_libphgeom_SOURCES = testexternals.cc
testexternals_libphgeom_LDADD   = libphgeom.la

testexternals_libphgeom_io_SOURCES = testexternals.cc
testexternals_libphgeom_io_LDADD   = libphgeom_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
