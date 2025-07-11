AUTOMAKE_OPTIONS = foreign

bin_PROGRAMS = \
  CaloCDB-GenStatus \
  CaloCDB-FilterDatasets

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
   `root-config --libs`

pkginclude_HEADERS = \
  genStatus.h \
  filter-datasets.h \
  geometry_constants.h \
  myUtils.h

lib_LTLIBRARIES = \
  libcalo_cdb.la

libcalo_cdb_la_SOURCES = \
  genStatus.cc \
  filter-datasets.cc \
  myUtils.cc

libcalo_cdb_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lcalo_io \
  -lcdbobjects \
  -lsphenixnpc \
  -lemcNoisyTowerFinder

CaloCDB_GenStatus_SOURCES = CaloCDB-GenStatus.cc
# CaloCDB_GenStatus_CXXFLAGS = -fsanitize=address
CaloCDB_GenStatus_LDADD = libcalo_cdb.la

CaloCDB_FilterDatasets_SOURCES = CaloCDB-FilterDatasets.cc
# CaloCDB_FilterDatasets_CXXFLAGS = -fsanitize=address
CaloCDB_FilterDatasets_LDADD = libcalo_cdb.la

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libcalo_cdb.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
