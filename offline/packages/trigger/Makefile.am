AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I${G4_MAIN}/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libcalotrigger_io.la \
   libcalotrigger.la

AM_CXXFLAGS = \
 -Werror -msse2

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs`

libcalotrigger_io_la_LIBADD = \
  -lphool

libcalotrigger_la_LIBADD = \
  -lfun4all \
  -lcalo_io -lcalo_util \
  -lCGAL \
  libcalotrigger_io.la

pkginclude_HEADERS = \
  CaloTriggerSim.h \
  CaloTriggerInfo.h \
  CaloTriggerInfo_v1.h

libcalotrigger_io_la_SOURCES = \
  CaloTriggerInfo_Dict.C \
  CaloTriggerInfo_v1.cc \
  CaloTriggerInfo_v1_Dict.C

libcalotrigger_la_SOURCES = \
  CaloTriggerSim.cc \
  CaloTriggerSim_Dict.C

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests
BUILT_SOURCES = testexternals.C

noinst_PROGRAMS = \
  testexternals_calotrigger_io \
  testexternals_calotrigger


testexternals_calotrigger_io_SOURCES = testexternals.C
testexternals_calotrigger_io_LDADD = libcalotrigger_io.la

testexternals_calotrigger_SOURCES = testexternals.C
testexternals_calotrigger_LDADD = libcalotrigger.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
