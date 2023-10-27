AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libg4mbd_io.la \
   libg4mbd.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4mbd_io_la_LIBADD = \
  -lphool

libg4mbd_la_LIBADD = \
  libg4mbd_io.la \
  -lmbd_io \
  -lg4detectors \
  -lSubsysReco

pkginclude_HEADERS = \
  MbdDigitization.h \
  MbdVertexFastSimReco.h


libg4mbd_io_la_SOURCES = \
  $(ROOTDICTS)

libg4mbd_la_SOURCES = \
  MbdDigitization.cc \
  MbdVertexFastSimReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternalsg4mbd_io \
  testexternalsg4mbd

BUILT_SOURCES = testexternals.cc

testexternalsg4mbd_io_SOURCES = testexternals.cc
testexternalsg4mbd_io_LDADD = libg4mbd_io.la

testexternalsg4mbd_SOURCES = testexternals.cc
testexternalsg4mbd_LDADD = libg4mbd_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
