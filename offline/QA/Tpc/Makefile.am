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
  Biquad.h \
  TpcChanQA.h \
  TpcFFT.h \
  TpcRawHitQA.h \
  TpcNoiseQA.h

lib_LTLIBRARIES = \
  libtpcqa.la

libtpcqa_la_SOURCES = \
  Biquad.cc \
  TpcChanQA.cc \
  TpcFFT.cc \
  TpcRawHitQA.cc \
  TpcNoiseQA.cc

libtpcqa_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lg4tpc \
  -ltrack_io \
  -ltrackbase_historic_io \
  -ltrack_reco \
  -lqautils

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libtpcqa.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
