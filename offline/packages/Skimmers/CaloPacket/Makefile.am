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
  CaloPacketSkimmer.h

lib_LTLIBRARIES = \
  libCaloPacketSkimmer.la

libCaloPacketSkimmer_la_SOURCES = \
  CaloPacketSkimmer.cc

libCaloPacketSkimmer_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lcalo_io \
  -lffarawobjects \
  -lg4jets \
  -lqautils


BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libCaloPacketSkimmer.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
