AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = libLiteCaloEvalTowSlope.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -lcalo_io \
  -lSubsysReco \
  -lfun4all \
  -lSubsysReco

libLiteCaloEvalTowSlope_la_SOURCES = \
  LiteCaloEval.cc \
  HCalCosmics.cc

pkginclude_HEADERS = \
  LiteCaloEval.h \
  HCalCosmics.h

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libLiteCaloEvalTowSlope.la

testexternals.cc:
	echo "this is a test.  $(ROOTSYS)/include"
	echo $(ROOTSYS)
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
