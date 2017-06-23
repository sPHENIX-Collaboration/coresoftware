AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir` \
  -I/gpfs/mnt/gpfs04/sphenix/user/nfeege/eicsw/HEPGen/install/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L/gpfs/mnt/gpfs04/sphenix/user/nfeege/eicsw/HEPGen/install/lib

AM_CXXFLAGS = -Wall -Werror

lib_LTLIBRARIES = libsHEPGen.la 

pkginclude_HEADERS = \
  sHEPGen.h

libsHEPGen_la_LDFLAGS = \
  ${AM_LDFLAGS} \
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
  sHEPGen.C \
  sHEPGen_Dict.C

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

sHEPGen_Dict.C: \
  sHEPGen.h \
  sHEPGenLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict*
	rm -f testexternals*
