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
  -lg4dst

pkginclude_HEADERS = \
  BaseTruthEval.h \
  CaloEvalStack.h \
  CaloTruthEval.h \
  CaloRawTowerEval.h \
  CaloRawClusterEval.h \
  CaloEvaluator.h \
  JetEvalStack.h \
  JetTruthEval.h \
  JetRecoEval.h \
  JetEvaluator.h \
  SvtxEvalStack.h \
  SvtxTruthEval.h \
  SvtxHitEval.h \
  SvtxClusterEval.h \
  SvtxTrackEval.h \
  SvtxVertexEval.h \
  SvtxEvaluator.h \
  MomentumEvaluator.h \
  PHG4DstCompressReco.h

#pkginclude_HEADERS = $(include_HEADERS)

libg4eval_la_SOURCES = \
  BaseTruthEval.C \
  CaloEvalStack.C \
  CaloEvalStack_Dict.C \
  CaloTruthEval.C \
  CaloTruthEval_Dict.C \
  CaloRawTowerEval.C \
  CaloRawTowerEval_Dict.C \
  CaloRawClusterEval.C \
  CaloRawClusterEval_Dict.C \
  CaloEvaluator.C \
  CaloEvaluator_Dict.C \
  JetEvalStack.C \
  JetEvalStack_Dict.C \
  JetTruthEval.C \
  JetTruthEval_Dict.C \
  JetRecoEval.C \
  JetRecoEval_Dict.C \
  JetEvaluator.C \
  JetEvaluator_Dict.C \
  PHG4DSTReader.cc \
  PHG4DSTReader_Dict.C \
  SvtxEvalStack.C \
  SvtxEvalStack_Dict.C \
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
  SvtxEvaluator_Dict.C \
  MomentumEvaluator.C \
  MomentumEvaluator_Dict.C \
  PHG4DstCompressReco.C \
  PHG4DstCompressReco_Dict.C

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
