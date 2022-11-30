##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libmicromegas_io.la \
  libmicromegas.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  CylinderGeomMicromegas.h \
  MicromegasClusterizer.h \
  MicromegasDefs.h \
  MicromegasTile.h

ROOTDICTS = \
  CylinderGeomMicromegas_Dict.cc \
  MicromegasTile_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  CylinderGeomMicromegas_Dict_rdict.pcm \
  MicromegasTile_Dict_rdict.pcm

# sources for io library
libmicromegas_io_la_SOURCES = \
  $(ROOTDICTS) \
  CylinderGeomMicromegas.cc \
  MicromegasDefs.cc

libmicromegas_io_la_LIBADD = \
  -lCLHEP \
  -lg4detectors_io \
  -lphool \
  -ltrack_io \
  -ltrackbase_historic_io

# sources for micromegas library
libmicromegas_la_SOURCES = \
  MicromegasClusterizer.cc

libmicromegas_la_LIBADD = \
  libmicromegas_io.la \
  -lphg4hit \
  -lSubsysReco

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_micromegas_io \
  testexternals_micromegas

testexternals_micromegas_io_SOURCES = testexternals.cc
testexternals_micromegas_io_LDADD = libmicromegas_io.la

testexternals_micromegas_SOURCES = testexternals.cc
testexternals_micromegas_LDADD = libmicromegas.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
