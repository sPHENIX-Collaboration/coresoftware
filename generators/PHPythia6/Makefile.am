AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -I$(includedir) -I$(OFFLINE_MAIN)/include -I`root-config --incdir`
AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib -L$(OPT_SPHENIX)/lib

AM_CXXFLAGS = -Wall -Werror

lib_LTLIBRARIES = libPHPythia6.la

include_HEADERS = \
  PHPythia6.h

libPHPythia6_la_LDFLAGS = ${AM_LDFLAGS} `root-config --libs`

libPHPythia6_la_LIBADD = \
  -lphool -lSubsysReco -lfun4all -lPythia6 -lEG -lEGPythia6 \
  -lphhepmc -lHepMC -lHepMCfio -lgsl -lgslcblas

libPHPythia6_la_SOURCES = \
  PHPythia6.C \
  PHPythia6_Dict.C

BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHPythia6.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

PHPythia6_Dict.C: \
  PHPythia6.h \
  PHPythia6LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict*
	rm -f testexternals*

