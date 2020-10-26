AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir` \
  -I$(OPT_SPHENIX)/HEPGenPlusPlus/include

lib_LTLIBRARIES = libsHEPGen.la

pkginclude_HEADERS = \
  sHEPGen.h

libsHEPGen_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/HEPGenPlusPlus/lib \
  `root-config --libs`

libsHEPGen_la_LIBADD = \
  -lgsl \
  -lgslcblas \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lboost_iostreams \
  -lVGG \
  -lhepgen \
  -lphhepmc \
  -lHepMC \
  -lCGAL

libsHEPGen_la_SOURCES = \
  sHEPGen.C

BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libsHEPGen.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
