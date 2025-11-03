AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  InttQa.h \
  InttRawHitQA.h \
  InttStreamQA.h \
  InttQaCommon.h \
  TH2INTT.h \
  InttMixupQA.h \
  bcocheck.h

lib_LTLIBRARIES = \
  libinttrawhitqa.la

libinttrawhitqa_la_SOURCES = \
  InttQa.cc \
  InttRawHitQA.cc \
  InttStreamQA.cc \
  TH2INTT.cc \
  InttMixupQA.cc \
  bcocheck.cc

libinttrawhitqa_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lg4tpc \
  -lg4intt \
  -lg4mvtx \
  -lmicromegas_io \
  -lmvtx_io \
  -lintt_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -ltrack_reco \
  -lqautils

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libinttrawhitqa.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
