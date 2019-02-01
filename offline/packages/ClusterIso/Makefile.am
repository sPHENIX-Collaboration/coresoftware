AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

pkginclude_HEADERS = \
  ClusterIso.h

libclusteriso_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

lib_LTLIBRARIES = \
  libclusteriso.la

if ! MAKEROOT6
ROOT5_DICTS = \
  ClusterIso_Dict.cc 
endif

libclusteriso_la_SOURCES = \
  ClusterIso.cc

libclusteriso_la_LIBADD = \
  -lcalo_io \
  -lg4vertex_io \
  -lSubsysReco

#  -ljetbackground
#  -lphool

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc


testexternals_LDADD = \
  libclusteriso.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f *Dict*  $(BUILT_SOURCES) *.pcm
