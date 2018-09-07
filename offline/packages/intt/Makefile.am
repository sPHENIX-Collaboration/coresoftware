AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libintt.la \
  libintt_io.la


AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include 

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

# sources for intt library
libintt_la_SOURCES = \
  InttDefs.cc \
  InttDefs_Dict.cc 

libintt_la_LIBADD = \
  -lfun4all \
  -ltrack_util \
  -ltrack_io

# sources for io library
libintt_io_la_SOURCES = \
  InttHit.cc \
  InttHit_Dict.cc

libintt_io_la_LIBADD = \
  libintt.la \
  -lfun4all \
  -ltrack_util \
  -ltrack_io 


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^


##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  InttDefs.h \
  InttHit.h



################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_intt_io \
  testexternals_intt

testexternals_intt_io_SOURCES = testexternals.cc
testexternals_intt_io_LDADD = libintt_io.la

testexternals_intt_SOURCES = testexternals.cc
testexternals_intt_LDADD = libintt.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* *.pcm
