AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

pkginclude_HEADERS = \
  TpcRS.h

lib_LTLIBRARIES = \
  libtpcresponse.la

libtpcresponse_la_SOURCES = \
  TpcRS.cc

libtpcresponse_la_LIBADD = \
  -lphg4hit \
  -lphool \
  -lSubsysReco \
  -ltpcrs \
  -lGeom \
  -lMathMore

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libtpcresponse.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
