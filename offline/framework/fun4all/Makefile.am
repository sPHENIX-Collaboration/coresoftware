AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  Fun4AllBase.h \
  Fun4AllDstInputManager.h \
  Fun4AllDstOutputManager.h \
  Fun4AllDummyInputManager.h \
  Fun4AllHistoBinDefs.h \
  Fun4AllHistoManager.h \
  Fun4AllInputManager.h \
  Fun4AllMemoryTracker.h \
  Fun4AllMonitoring.h \
  Fun4AllNoSyncDstInputManager.h \
  Fun4AllOutputManager.h \
  Fun4AllReturnCodes.h \
  Fun4AllServer.h \
  Fun4AllSyncManager.h \
  Fun4AllUtils.h \
  PHTFileServer.h \
  SubsysReco.h \
  TDirectoryHelper.h

lib_LTLIBRARIES = \
  libSubsysReco.la \
  libTDirectoryHelper.la \
  libfun4all.la

libTDirectoryHelper_la_SOURCES = \
  TDirectoryHelper.cc

libTDirectoryHelper_la_LDFLAGS = \
    `root-config --libs`

libfun4all_la_SOURCES = \
  Fun4AllDstInputManager.cc \
  Fun4AllDstOutputManager.cc \
  Fun4AllDummyInputManager.cc \
  Fun4AllHistoManager.cc \
  Fun4AllInputManager.cc \
  Fun4AllMonitoring.cc \
  Fun4AllMemoryTracker.cc \
  Fun4AllNoSyncDstInputManager.cc \
  Fun4AllOutputManager.cc \
  Fun4AllServer.cc \
  Fun4AllSyncManager.cc \
  Fun4AllUtils.cc \
  PHTFileServer.cc

libfun4all_la_LIBADD = \
  libSubsysReco.la \
  libTDirectoryHelper.la \
  -lboost_filesystem \
  -lFROG \
  -lffaobjects \
  -lphool

libSubsysReco_la_SOURCES = \
  Fun4AllBase.cc

bin_SCRIPTS = \
  CreateSubsysRecoModule.pl

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_fun4all \
  testexternals_subsysreco \
  testexternals_tdirectoryhelper

testexternals_fun4all_SOURCES = testexternals.cc
testexternals_fun4all_LDADD   = libfun4all.la

testexternals_subsysreco_SOURCES = testexternals.cc
testexternals_subsysreco_LDADD   = libSubsysReco.la

testexternals_tdirectoryhelper_SOURCES = testexternals.cc
testexternals_tdirectoryhelper_LDADD = libTDirectoryHelper.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
