AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(OPT_SPHENIX)/include \
  -isystem$(ROOTSYS)/include 

lib_LTLIBRARIES = \
  libfun4cal.la

libfun4cal_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  -lfun4all \
  -lodbc++ \
  -lpdbcalBase \
  -lphool

pkginclude_HEADERS = \
  Fun4CalDBCodes.h \
  Fun4CalHistoBinDefs.h \
  CalReco.h \
  Fun4CalServer.h

libfun4cal_la_SOURCES = \
  CalReco.cc \
  Fun4CalServer.cc

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libfun4cal.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@
