AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I$(OFFLINE_MAIN)/include \
  -I${G4_MAIN}/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libg4eval.la

libg4eval_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4eval_la_LIBADD = \
  -lcalo_io \
  -lCLHEP \
  -lfun4all \
  -lg4detectors_io \
  -lg4jets_io \
  -lg4vertex_io \
  -ltrackbase_historic_io \
  -ltrack_io \
  -ltpc_io \
  -lmvtx_io \
  -lintt_io \
  -lmicromegas_io \
  -lphg4hit \
  -lgsl \
  -lgslcblas

pkginclude_HEADERS = \
  BaseTruthEval.h \
  CaloEvalStack.h \
  CaloEvaluator.h \
  CaloRawClusterEval.h \
  CaloRawTowerEval.h \
  CaloTruthEval.h \
  EventCounter_hp.h \
  JetEvalStack.h \
  JetEvaluator.h \
  JetRecoEval.h \
  JetTruthEval.h \
  MicromegasEvaluator_hp.h \
  MomentumEvaluator.h \
  PHG4DstCompressReco.h \
  PHG4DSTReader.h \
  SimEvaluator_hp.h \
  SvtxEvalStack.h \
  SvtxEvaluator.h \
  SvtxHitEval.h \
  SvtxClusterEval.h \
  SvtxTrackEval.h \
  SvtxTruthEval.h \
  SvtxVertexEval.h \
  TrackingEvaluator_hp.h

if ! MAKEROOT6
  ROOT5DICTS = \
  CaloEvalStack_Dict.cc \
  CaloTruthEval_Dict.cc \
  CaloRawTowerEval_Dict.cc \
  CaloRawClusterEval_Dict.cc \
  CaloEvaluator_Dict.cc \
  JetEvalStack_Dict.cc \
  JetTruthEval_Dict.cc \
  JetRecoEval_Dict.cc \
  JetEvaluator_Dict.cc \
  PHG4DSTReader_Dict.cc \
  SvtxEvalStack_Dict.cc \
  SvtxTruthEval_Dict.cc \
  SvtxHitEval_Dict.cc \
  SvtxClusterEval_Dict.cc \
  SvtxTrackEval_Dict.cc \
  SvtxVertexEval_Dict.cc \
  SvtxEvaluator_Dict.cc \
  MomentumEvaluator_Dict.cc \
  PHG4DstCompressReco_Dict.cc
endif

libg4eval_la_SOURCES = \
  $(ROOT5DICTS) \
  BaseTruthEval.cc \
  CaloEvalStack.cc \
  CaloTruthEval.cc \
  CaloRawTowerEval.cc \
  CaloRawClusterEval.cc \
  CaloEvaluator.cc \
  EventCounter_hp.cc \
  JetEvalStack.cc \
  JetTruthEval.cc \
  JetRecoEval.cc \
  JetEvaluator.cc \
  MicromegasEvaluator_hp_Dict.cc \
  MicromegasEvaluator_hp.cc \
  PHG4DSTReader.cc \
  SimEvaluator_hp.cc \
  SimEvaluator_hp_Dict.cc \
  SvtxEvalStack.cc \
  SvtxTruthEval.cc \
  SvtxHitEval.cc \
  SvtxClusterEval.cc \
  SvtxTrackEval.cc \
  SvtxVertexEval.cc \
  SvtxEvaluator.cc \
  TrackingEvaluator_hp.cc \
  TrackingEvaluator_hp_Dict.cc \
  MomentumEvaluator.cc \
  PHG4DstCompressReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libg4eval.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
