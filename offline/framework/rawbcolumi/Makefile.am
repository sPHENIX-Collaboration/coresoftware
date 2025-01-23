AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  SingleGl1PoolInputv2.h \
  SingleStreamingInputv2.h \
  Fun4AllStreamingLumiCountingInputManager.h

lib_LTLIBRARIES = \
  librawbcolumi.la

# source for mvtx decoder library
librawbcolumi_la_SOURCES = \
  SingleGl1PoolInputv2.cc \
  SingleStreamingInputv2.cc \
  Fun4AllStreamingLumiCountingInputManager.cc

librawbcolumi_la_LIBADD = \
  -lffarawobjects \
  -lfun4all \
  -lEvent \
  -lphoolraw \
  -lqautils

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = librawbcolumi.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
