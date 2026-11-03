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
  BcoRangeCheck.h \
  ClockDiffCheck.h \
  DumpPacket.h \
  EventCombiner.h \
  EventNumberCheck.h \
  Gl1BcoDump.h \
  Gl1Check.h \
  InttCheck.h \
  InttBcoDump.h \
  InttGl1Check.h \
  LL1Check.h \
  StreamingCheck.h \
  TpcBcoDump.h \
  TpcCheck.h

lib_LTLIBRARIES = \
  libffarawmodules.la

libffarawmodules_la_SOURCES = \
  BcoRangeCheck.cc \
  ClockDiffCheck.cc \
  DumpPacket.cc \
  EventCombiner.cc \
  EventNumberCheck.cc \
  Gl1BcoDump.cc \
  Gl1Check.cc \
  InttBcoDump.cc \
  InttCheck.cc \
  InttGl1Check.cc \
  LL1Check.cc \
  StreamingCheck.cc \
  TpcBcoDump.cc \
  TpcCheck.cc

libffarawmodules_la_LIBADD = \
  -lEvent \
  -lffarawobjects \
  -lfun4allraw \
  -lphool \
  -lSubsysReco

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libffarawmodules.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
