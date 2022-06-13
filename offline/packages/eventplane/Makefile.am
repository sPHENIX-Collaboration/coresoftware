AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  EpFinder.h \
  EpInfo.h \
  EpFinderReco.h

lib_LTLIBRARIES = \
  libeventplane.la

libeventplane_la_SOURCES = \
  EpFinder.cc \
  EpInfo.cc \
  EpFinderReco.cc

libeventplane_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lcentrality_io \
  -lphg4hit \
  -lcalo_io

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libeventplane.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
