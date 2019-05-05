AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libSeamstress.la

pkginclude_HEADERS = \
  Seamstress.h \
  Needle.h \
  Pincushion.h

libSeamstress_la_SOURCES = \
  Seamstress.cpp

libSeamstress_la_LIBADD = \
  -lpthread

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libSeamstress.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
