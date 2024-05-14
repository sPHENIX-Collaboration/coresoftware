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
  StructureinJets.h \
  TrksInJetQA.h \
  TrksInJetQAHist.h \
  TrksInJetQATypes.h \
  TrksInJetQAConfig.h \
  TrksInJetQABaseFiller.h \
  TrksInJetQAInJetFiller.h \
  TrksInJetQAInclusiveFiller.h \
  TrksInJetQABaseManager.h \
  TrksInJetQAHitManager.h \
  TrksInJetQAClustManager.h \
  TrksInJetQATrkManager.h \
  TrksInJetQAJetManager.h \
  JetSeedCount.h

libjetqa_la_SOURCES = \
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
  JetSeedCount.cc

libjetqa_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lcalo_io \
  -lcentrality_io \
  -lparticleflow_io \
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

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f $(BUILT_SOURCES)
