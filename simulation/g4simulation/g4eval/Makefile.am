AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I$(OFFLINE_MAIN)/include \
  -I$(G4_MAIN)/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libg4eval.la

AM_CXXFLAGS = -Werror -msse2

libg4eval_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4eval_la_LIBADD = \
  -lfun4all \
  -ltrackbase_historic_io \
  -lg4dst 

pkginclude_HEADERS = \
  BaseTruthEval.h \
  CaloEvalStack.h \
  CaloEvaluator.h \
  CaloRawClusterEval.h \
  CaloRawTowerEval.h \
  CaloTruthEval.h \
  JetEvalStack.h \
  JetEvaluator.h \
  JetRecoEval.h \
  JetTruthEval.h \
  MomentumEvaluator.h \
  PHG4DstCompressReco.h \
  PHG4DSTReader.h \
  SvtxEvalStack.h \
  SvtxEvaluator.h \
  SvtxHitEval.h \
  SvtxClusterEval.h \
  SvtxTrackEval.h \
  SvtxTruthEval.h \
  SvtxVertexEval.h

if ! MAKEROOT6
  ROOT5DICTS = \
  CaloEvalStack_Dict.C \
  CaloTruthEval_Dict.C \
  CaloRawTowerEval_Dict.C \
  CaloRawClusterEval_Dict.C \
  CaloEvaluator_Dict.C \
  JetEvalStack_Dict.C \
  JetTruthEval_Dict.C \
  JetRecoEval_Dict.C \
  JetEvaluator_Dict.C \
  PHG4DSTReader_Dict.C \
  SvtxEvalStack_Dict.C \
  SvtxTruthEval_Dict.C \
  SvtxHitEval_Dict.C \
  SvtxClusterEval_Dict.C \
  SvtxTrackEval_Dict.C \
  SvtxVertexEval_Dict.C \
  SvtxEvaluator_Dict.C \
  MomentumEvaluator_Dict.C \
  PHG4DstCompressReco_Dict.C
endif

libg4eval_la_SOURCES = \
  BaseTruthEval.C \
  CaloEvalStack.C \
  CaloTruthEval.C \
  CaloRawTowerEval.C \
  CaloRawClusterEval.C \
  CaloEvaluator.C \
  JetEvalStack.C \
  JetTruthEval.C \
  JetRecoEval.C \
  JetEvaluator.C \
  PHG4DSTReader.cc \
  SvtxEvalStack.C \
  SvtxTruthEval.C \
  SvtxHitEval.C \
  SvtxClusterEval.C \
  SvtxTrackEval.C \
  SvtxVertexEval.C \
  SvtxEvaluator.C \
  MomentumEvaluator.C \
  PHG4DstCompressReco.C 

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

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
	rm -f *Dict* testexternals.C *.pcm
