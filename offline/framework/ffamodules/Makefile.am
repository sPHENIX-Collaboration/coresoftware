AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libffamodules.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libffamodules_la_LIBADD = \
  -lfun4all \
  -lffaobjects \
  -lSubsysReco


pkginclude_HEADERS = \
  SyncReco.h

libffamodules_la_SOURCES = \
  SyncReco.cc

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libffamodules.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
