AUTOMAKE_OPTIONS = foreign

INCLUDES = \
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
  -lcemc_io \
  -lg4jets_io \
  -lg4eval

pkginclude_HEADERS = 

#pkginclude_HEADERS = $(include_HEADERS)

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
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

%_Dict.cpp: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

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
