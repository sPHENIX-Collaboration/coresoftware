AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libfun4allutils.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

libfun4allutils_la_LIBADD = \
  -lcdbobjects \
  -lfun4all \
  -lSubsysReco

pkginclude_HEADERS = \
  TimerStats.h

libfun4allutils_la_SOURCES = \
  TimerStats.cc

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libfun4allutils.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
