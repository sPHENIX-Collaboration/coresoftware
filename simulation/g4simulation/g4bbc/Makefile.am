AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libg4bbc_io.la \
   libg4bbc.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4bbc_io_la_LIBADD = \
  -lphool

libg4bbc_la_LIBADD = \
  libg4bbc_io.la \
  -lbbc_io \
  -lg4detectors \
  -lSubsysReco

pkginclude_HEADERS = \
  BbcDigitization.h \
  BbcVertexFastSimReco.h


libg4bbc_io_la_SOURCES = \
  $(ROOTDICTS)

libg4bbc_la_SOURCES = \
  BbcDigitization.cc \
  BbcVertexFastSimReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternalsg4bbc_io \
  testexternalsg4bbc

BUILT_SOURCES = testexternals.cc

testexternalsg4bbc_io_SOURCES = testexternals.cc
testexternalsg4bbc_io_LDADD = libg4bbc_io.la

testexternalsg4bbc_SOURCES = testexternals.cc
testexternalsg4bbc_LDADD = libg4bbc_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
