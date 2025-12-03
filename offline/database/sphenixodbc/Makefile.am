AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(OPT_SPHENIX)/include

lib_LTLIBRARIES = libsphenixodbc.la

bin_PROGRAMS = testODBC

libsphenixodbc_la_LIBADD = \
  -L$(OPT_SPHENIX)/lib \
  -lodbc++

pkginclude_HEADERS = \
  ODBCInterface.h

libsphenixodbc_la_SOURCES = \
  ODBCInterface.cc

noinst_PROGRAMS = testexternals

testODBC_SOURCES = testODBC.cc
testODBC_LDADD = libsphenixodbc.la

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = testexternals.cc

testexternals_LDADD = \
  libsphenixodbc.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
