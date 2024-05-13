AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libuspin.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include

pkginclude_HEADERS = \
  SpinDBInput.h \
  SpinDBContent.h \
  SpinDBOutput.h

libuspin_la_SOURCES = \
  SpinDBInput.cc \
  SpinDBContent.cc \
  SpinDBOutput.cc \
  SpinDBUtil_Dict.C

libuspin_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  -L$(ROOTSYS)/lib \
  -lCore \
  -lodbc++

noinst_HEADERS =  \
  SpinDBUtilLinkDef.h

pcmdir = $(libdir)

nobase_dist_pcm_DATA = \
  SpinDBUtil_Dict_rdict.pcm

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

SpinDBUtil_Dict.C: SpinDBContent.h SpinDBInput.h SpinDBOutput.h SpinDBUtilLinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
