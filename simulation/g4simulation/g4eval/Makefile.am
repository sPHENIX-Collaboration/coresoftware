AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I$(OFFLINE_MAIN)/include \
  -I$(G4_MAIN)/include \
  -I$(G4_MAIN)/include/Geant4 \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libg4eval_io.la \
   libg4eval.la

AM_CXXFLAGS = -Wall -Werror -msse2

libg4eval_io_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool

libg4eval_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs` \
  -lfun4all \
  -lg4hough \
  -lHelixHough \
  -lFitNewton \
  -lcemc \
  libg4eval_io.la

pkginclude_HEADERS = \
  PHG4Evaluator.h \
  PHG4CalEvaluator.h \
  PHG4DSTReader.h \
  EvalLinks.h \
  EvalLinksV1.h \
  SvtxClusterEval.h

#pkginclude_HEADERS = $(include_HEADERS)

libg4eval_io_la_SOURCES = \
  EvalLinks.C \
  EvalLinks_Dict.C \
  EvalLinksV1.C \
  EvalLinksV1_Dict.C

libg4eval_la_SOURCES = \
  PHG4Evaluator.C \
  PHG4Evaluator_Dict.C \
  PHG4CalEvaluator.C \
  PHG4CalEvaluator_Dict.C \
  PHG4DSTReader.cc \
  PHG4DSTReader_Dict.C \
  SvtxClusterEval.C \
  SvtxClusterEval_Dict.C

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
  libg4eval_io.la \
  libg4eval.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
