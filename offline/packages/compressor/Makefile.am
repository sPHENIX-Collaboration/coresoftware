AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libcompressor.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include 

libcompressor_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

pkginclude_HEADERS = \
  compressor.h

libcompressor_la_SOURCES = \
  compress_clu_res_float32.cc

################################################
# linking test to make sure we do not have unresolved symbols

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libcompressor.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f testexternals.cc
