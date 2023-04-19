AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libffamodules.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

libffamodules_la_LIBADD = \
  -lfun4all \
  -lffaobjects \
  -lphhepmc_io \
  -lsphenixnpc \
  -lSubsysReco

pkginclude_HEADERS = \
  CDBInterface.h \
  FlagHandler.h \
  HeadReco.h \
  SyncReco.h

libffamodules_la_SOURCES = \
  CDBInterface.cc \
  FlagHandler.cc \
  HeadReco.cc \
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
