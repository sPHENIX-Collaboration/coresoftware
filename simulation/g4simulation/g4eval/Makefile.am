AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I$(OFFLINE_MAIN)/include \
  -I$(G4_MAIN)/include \
  -I$(G4_MAIN)/include/Geant4 \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libg4eval.la

AM_CXXFLAGS = -Wall -Werror -msse2

libg4eval_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs` \
  -lfun4all \
  -lg4hough \
  -lHelixHough \
  -lFitNewton \
  -lcemc 

pkginclude_HEADERS = \
  CaloTruthEval.h \
  CaloRawTowerEval.h \
  PHG4Evaluator.h \
  PHG4CalEvaluator.h \
  PHG4DSTReader.h \
  SvtxTruthEval.h \
  SvtxHitEval.h \
  SvtxClusterEval.h \
  SvtxTrackEval.h \
  SvtxVertexEval.h \
  SvtxEvaluator.h

#pkginclude_HEADERS = $(include_HEADERS)

libg4eval_la_SOURCES = \
  CaloTruthEval.C \
  CaloTruthEval_Dict.C \
  CaloRawTowerEval.C \
  CaloRawTowerEval_Dict.C \
  PHG4Evaluator.C \
  PHG4Evaluator_Dict.C \
  PHG4CalEvaluator.C \
  PHG4CalEvaluator_Dict.C \
  PHG4DSTReader.cc \
  PHG4DSTReader_Dict.C \
  SvtxTruthEval.C \
  SvtxTruthEval_Dict.C \
  SvtxHitEval.C \
  SvtxHitEval_Dict.C \
  SvtxClusterEval.C \
  SvtxClusterEval_Dict.C \
  SvtxTrackEval.C \
  SvtxTrackEval_Dict.C \
  SvtxVertexEval.C \
  SvtxVertexEval_Dict.C \
  SvtxEvaluator.C \
  SvtxEvaluator_Dict.C

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
  libg4eval.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
