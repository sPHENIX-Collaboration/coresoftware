AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libqa_modules.la

AM_CXXFLAGS = -Wall -Werror -msse2

libqa_modules_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libqa_modules_la_LIBADD = \
  -lfun4all \
  -lphg4hit \
  -lg4detectors_io \
  -lg4hough_io \
  -lg4hough \
  -lcalo_io -lcalo_util \
  -lg4jets_io \
  -lg4eval

libqa_modules_la_SOURCES = \
  QAHistManagerDef.C \
  QAHistManagerDef_Dict.C \
  QAG4SimulationCalorimeter.C \
  QAG4SimulationCalorimeter_Dict.C \
  QAG4SimulationCalorimeterSum.C \
  QAG4SimulationCalorimeterSum_Dict.C \
  QAG4SimulationJet.C \
  QAG4SimulationJet_Dict.C 

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libqa_modules.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* testexternals.C
