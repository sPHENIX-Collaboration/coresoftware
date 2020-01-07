AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libphoolraw.la


AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libphoolraw_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

libphoolraw_la_LIBADD = \
  -lphool \
  -lEvent

ROOTDICTS =

if MAKEROOT6
pcmdir = $(libdir)
nobase_dist_pcm_DATA =

else
ROOT5DICTS =
endif

libphoolraw_la_SOURCES = \
  PHRawDataNode.cc \
  PHRawOManager.cc

pkginclude_HEADERS =  \
  PHRawOManager.h

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals


testexternals_SOURCES = testexternals.cc
testexternals_LDADD = \
  libphoolraw.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
