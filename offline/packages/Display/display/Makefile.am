
AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libpheve_display.la

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib `root-config --evelibs`

AM_CXXFLAGS = `root-config --cflags`

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include -I$(OPT_PHENIX)/include

include_HEADERS = \
  PHEveDisplay.h

libpheve_display_la_LIBADD = -lboost_thread 

#-lmfm

libpheve_display_la_SOURCES = \
  PHEveDisplay.cxx

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libpheve_display.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)

%Dict.C : %.h %LinkDef.h
	rootcint -f $@ -c -D__linux__ $(DEFAULT_INCLUDES) $(INCLUDES) $^

