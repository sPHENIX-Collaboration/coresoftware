AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I${G4_MAIN}/include \
  -I${G4_MAIN}/include/Geant4 \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libcalotrigger_io.la \
   libcalotrigger.la

AM_CXXFLAGS = \
 -Wall -Werror -msse2

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs`

libcalotrigger_io_la_LIBADD = \
  -lphool

libcalotrigger_la_LIBADD = \
  -lfun4all \
  -lcemc_io \
  -lCGAL \
  libcalotrigger_io.la

pkginclude_HEADERS = \
  CaloTriggerSim.h \
  CaloTriggerInfo.h

#pkginclude_HEADERS = $(include_HEADERS)

libcalotrigger_io_la_SOURCES = \
  CaloTriggerInfo.C \
  CaloTriggerInfo_Dict.C

libcalotrigger_la_SOURCES = \
  CaloTriggerSim.C \
  CaloTriggerSim_Dict.C

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

%_Dict.cpp: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libcalotrigger_io.la \
  libcalotrigger.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* testexternals.C
