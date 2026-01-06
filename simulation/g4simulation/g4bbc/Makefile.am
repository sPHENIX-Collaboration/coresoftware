AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libg4mbd.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4mbd_la_LIBADD = \
  -lmbd_io \
  -lg4detectors \
  -lSubsysReco \
  -lglobalvertex_io

pkginclude_HEADERS = \
  MbdDigitization.h \
  MbdVertexFastSimReco.h

libg4mbd_la_SOURCES = \
  MbdDigitization.cc \
  MbdVertexFastSimReco.cc

################################################
# linking tests

noinst_PROGRAMS = \
  testexternalsg4mbd

BUILT_SOURCES = testexternals.cc

testexternalsg4mbd_SOURCES = testexternals.cc
testexternalsg4mbd_LDADD = libg4mbd.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
