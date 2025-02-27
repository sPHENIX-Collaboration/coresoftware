AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir` \
  -isystem$(OPT_SPHENIX)/include

lib_LTLIBRARIES = \
   libcalotrigger_io.la \
   libcalotrigger.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libcalotrigger_io_la_LIBADD = \
  -lphool \
  -lodbc++ \
  -lcalo_io \
  -lmbd_io \
  -lzdcinfo_io

libcalotrigger_la_LIBADD = \
  libcalotrigger_io.la \
  -lmbd_io \
  -lzdcinfo_io \
  -lodbc++ \
  -lglobalvertex \
  -lcalo_reco \
  -lffamodules \
  -lSubsysReco \
  -lphool 

pkginclude_HEADERS = \
  MinimumBiasClassifier.h \
  MinimumBiasInfo.h \
  MinimumBiasInfov1.h \
  CaloTriggerSim.h \
  CaloTriggerInfo.h \
  CaloTriggerInfov1.h \
  LL1Defs.h \
  LL1Out.h \
  LL1Outv1.h \
  TriggerAnalyzer.h \
  CaloTriggerEmulator.h \
  MBDTriggerEmulator.h \
  LL1PacketGetter.h \
  TriggerDefs.h \
  TriggerRunInfoReco.h \
  TriggerRunInfo.h \
  TriggerRunInfov1.h \
  TriggerPrimitive.h \
  TriggerPrimitivev1.h \
  TriggerPrimitiveContainer.h \
  TriggerPrimitiveContainerv1.h


ROOTDICTS = \
  LL1Out_Dict.cc \
  LL1Outv1_Dict.cc \
  TriggerRunInfo_Dict.cc \
  TriggerRunInfov1_Dict.cc \
  TriggerPrimitive_Dict.cc \
  TriggerPrimitivev1_Dict.cc \
  TriggerPrimitiveContainer_Dict.cc \
  TriggerPrimitiveContainerv1_Dict.cc \
  CaloTriggerInfo_Dict.cc \
  CaloTriggerInfov1_Dict.cc \
  MinimumBiasInfo_Dict.cc \
  MinimumBiasInfov1_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libcalotrigger_io_la_SOURCES = \
  $(ROOTDICTS) \
  TriggerRunInfo.cc \
  TriggerRunInfov1.cc \
  LL1Out.cc \
  LL1Outv1.cc \
  TriggerPrimitive.cc \
  TriggerPrimitivev1.cc \
  TriggerPrimitiveContainer.cc \
  TriggerPrimitiveContainerv1.cc \
  TriggerDefs.cc \
  CaloTriggerInfov1.cc \
  MinimumBiasInfov1.cc

libcalotrigger_la_SOURCES = \
  TriggerRunInfoReco.cc \
  LL1PacketGetter.cc \
  TriggerAnalyzer.cc\
  CaloTriggerEmulator.cc \
  MBDTriggerEmulator.cc \
  CaloTriggerSim.cc \
  MinimumBiasClassifier.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_calotrigger_io \
  testexternals_calotrigger

testexternals_calotrigger_io_SOURCES = testexternals.cc
testexternals_calotrigger_io_LDADD = libcalotrigger_io.la

testexternals_calotrigger_SOURCES = testexternals.cc
testexternals_calotrigger_LDADD = libcalotrigger.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
