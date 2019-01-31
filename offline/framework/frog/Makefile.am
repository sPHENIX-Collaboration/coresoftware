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
  FROG.h

if ! MAKEROOT6
ROOT5DICTS = \
  FROG_Dict.cc
endif

libFROG_la_SOURCES = \
  $(ROOT5DICTS) \
  FROG.cc

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = testexternals.cc

testexternals_LDADD = \
  libFROG.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
