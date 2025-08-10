AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libuspin.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  `root-config --libs`

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include

pkginclude_HEADERS = \
  SpinDBInput.h \
  SpinDBContent.h \
  SpinDBContentv1.h \
  SpinDBOutput.h \
  SpinDBNode.h

libuspin_la_SOURCES = \
  $(ROOTDICTS) \
  SpinDBInput.cc \
  SpinDBContent.cc \
  SpinDBContentv1.cc \
  SpinDBOutput.cc \
  SpinDBNode.cc

libuspin_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  -lodbc++ \
  -lphool \
  -lfun4all \
  -lffaobjects \
  -lffarawobjects

ROOTDICTS = \
  SpinDBContent_Dict.cc \
  SpinDBContentv1_Dict.cc 

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libuspin.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
