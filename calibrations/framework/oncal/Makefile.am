AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OPT_SPHENIX)/include \
  -isystem$(ROOTSYS)/include 

lib_LTLIBRARIES = \
  liboncal.la

liboncal_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  -lfun4all \
  -lodbc++ \
  -lpdbcalBase \
  -lphool

pkginclude_HEADERS = \
  OnCalDBCodes.h \
  OnCalHistoBinDefs.h \
  OnCal.h \
  OnCalServer.h

liboncal_la_SOURCES = \
  OnCal.cc \
  OnCalServer.cc

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  liboncal.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@
