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

# Headers installed for use in macros and by other packages
pkginclude_HEADERS = \
  GenStatus.h \
  geometry_constants.h

# Headers used only for building this library and its binaries
noinst_HEADERS = \
  FilterDatasets.h \
  myUtils.h

lib_LTLIBRARIES = \
  libcalo_cdb.la

libcalo_cdb_la_SOURCES = \
  GenStatus.cc \
  FilterDatasets.cc \
  myUtils.cc

libcalo_cdb_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lcalo_io \
  -lcdbobjects \
  -lsphenixnpc \
  -lffamodules \
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
