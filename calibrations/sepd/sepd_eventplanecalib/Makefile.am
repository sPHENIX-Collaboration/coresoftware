AUTOMAKE_OPTIONS = foreign

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
  sEPD_TreeGen.h \
  QVecCalib.h \
  QVecDefs.h

lib_LTLIBRARIES = \
  libsepd_eventplanecalib.la

ROOTDICTS = \
  EventPlaneData_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

# EventPlaneData is a locally used root i/o object - no need to create an io library
libsepd_eventplanecalib_la_SOURCES = \
  $(ROOTDICTS) \
  EventPlaneData.cc \
  sEPD_TreeGen.cc \
  QVecCalib.cc

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

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

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
