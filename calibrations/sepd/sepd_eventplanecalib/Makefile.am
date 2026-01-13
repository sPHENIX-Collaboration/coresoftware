AUTOMAKE_OPTIONS = foreign

bin_PROGRAMS = \
  GenQVecCalib \
  GenQVecCDB

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
   `root-config --libs`

pkginclude_HEADERS = \
  sEPD_TreeGen.h \
  QVecCalib.h \
  QVecCDB.h \
  QVecDefs.h

lib_LTLIBRARIES = \
  libsepd_eventplanecalib.la

libsepd_eventplanecalib_la_SOURCES = \
  sEPD_TreeGen.cc \
  QVecCalib.cc \
  QVecCDB.cc

libsepd_eventplanecalib_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lcentrality_io \
  -lfun4all \
  -lffamodules \
  -lglobalvertex_io \
  -lcalotrigger_io \
  -lcalotrigger \
  -lcdbobjects \
  -lepd_io

GenQVecCalib_SOURCES = GenQVecCalib.cc
# GenQVecCalib_CXXFLAGS = -fsanitize=address
GenQVecCalib_LDADD = libsepd_eventplanecalib.la

GenQVecCDB_SOURCES = GenQVecCDB.cc
# GenQVecCDB_CXXFLAGS = -fsanitize=address
GenQVecCDB_LDADD = libsepd_eventplanecalib.la

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libsepd_eventplanecalib.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
