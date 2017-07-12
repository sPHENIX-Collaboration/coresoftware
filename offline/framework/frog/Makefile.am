AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OPT_SPHENIX)/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = libFROG.la

libFROG_la_LDFLAGS = \
  `root-config --libs`

libFROG_la_LIBADD = \
  -L$(libdir) \
  -L$(OPT_SPHENIX)/lib \
  -lodbc++

noinst_HEADERS = FROGLinkDef.h

pkginclude_HEADERS =    \
  pgsearch.h \
  dCachesearch.h  \
  FROG.h

libFROG_la_SOURCES = \
  pgsearch.cc  dCachesearch.cc \
  FROG.cc       \
  FROG_Dict.cc

BUILT_SOURCES = testexternals.cc
noinst_PROGRAMS = testexternals
testexternals_SOURCES = testexternals.cc

testexternals_LDADD = \
  libFROG.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating CINT dictionaries from class headers.
FROG_Dict.cc:   \
  pgsearch.h   \
  FROG.h       \
  FROGLinkDef.h
	rootcint -f $@ @CINTDEFS@ -c  $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict* *.pcm
