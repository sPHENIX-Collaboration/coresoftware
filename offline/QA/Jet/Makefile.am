AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libjetqa.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

pkginclude_HEADERS = \
  CaloStatusMapper.h \
  CaloStatusMapperDefs.h \
  ConstituentsinJets.h \
  DijetQA.h \
  EMClusterKinematics.h \
  JetKinematicCheck.h \
  JetQADefs.h \
  JetSeedCount.h \
  PhotonJetsKinematics.h \
  RhosinEvent.h \
  StructureinJets.h \
  TrksInJetQA.h \
  TrksInJetQAHist.h \
  TrksInJetQADefs.h \
  TrksInJetQAConfig.h \
  TrksInJetQABaseFiller.h \
  TrksInJetQAInJetFiller.h \
  TrksInJetQAInclusiveFiller.h \
  TrksInJetQABaseManager.h \
  TrksInJetQAHitManager.h \
  TrksInJetQAClustManager.h \
  TrksInJetQATrkManager.h \
  TrksInJetQAJetManager.h \
  UEvsEtaCentrality.h

libjetqa_la_SOURCES = \
  CaloStatusMapper.cc \
  ConstituentsinJets.cc \
  DijetQA.cc \
  EMClusterKinematics.cc \
  JetKinematicCheck.cc \
  JetSeedCount.cc \
  PhotonJetsKinematics.cc \
  RhosinEvent.cc \
  StructureinJets.cc \
  TrksInJetQA.cc \
  TrksInJetQABaseFiller.cc \
  TrksInJetQAInJetFiller.cc \
  TrksInJetQAInclusiveFiller.cc \
  TrksInJetQABaseManager.cc \
  TrksInJetQAHitManager.cc \
  TrksInJetQAClustManager.cc \
  TrksInJetQATrkManager.cc \
  TrksInJetQAJetManager.cc \
  UEvsEtaCentrality.cc

libjetqa_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lcalo_io \
  -lcalotrigger \
  -ljetbackground\
  -lg4jets \
  -lcentrality_io \
  -lparticleflow_io \
  -lffarawobjects \
  -lfun4all \
  -lphool \
  -lg4eval \
  -lparticleflow_io \
  -lSubsysReco \
  -lqautils

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libjetqa.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)


