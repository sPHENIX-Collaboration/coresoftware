AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libphgeom.la

pkginclude_HEADERS = \
	PHGeomTGeo.h

noinst_HEADERS = \
	PHGeomManager.h

libphgeom_la_SOURCES = \
  PHGeomTGeo.cc \
  PHGeomTGeo_Dict.cc \
  PHGeomManager.cc \
  PHGeomManager_Dict.cc 
 
libphgeom_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool
  
noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libphgeom.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^


clean-local:
	rm -f *Dict*
