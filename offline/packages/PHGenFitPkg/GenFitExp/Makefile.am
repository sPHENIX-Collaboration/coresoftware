AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

lib_LTLIBRARIES = \
  libgenfit2exp.la

pkginclude_HEADERS = \
  Field.h

libgenfit2exp_la_SOURCES = \
  Field.cc

libgenfit2exp_la_LIBADD = \
  -lgenfit2 \
  -lphfield

BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libgenfit2exp.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

testexternals_SOURCES = testexternals.C
