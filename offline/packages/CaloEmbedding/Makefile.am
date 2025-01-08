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
  caloTowerEmbed.h \
  CopyIODataNodes.h \
  HepMCCollisionVertex.h

lib_LTLIBRARIES = \
  libCaloEmbedding.la

libCaloEmbedding_la_SOURCES = \
  caloTowerEmbed.cc \
  CopyIODataNodes.cc \
  HepMCCollisionVertex.cc

libCaloEmbedding_la_LIBADD = \
  -lcalo_io \
  -lcalotrigger_io \
  -lcentrality_io \
  -lfun4all \
  -lglobalvertex_io \
  -lphhepmc_io \
  -lphool \
  -lSubsysReco

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libCaloEmbedding.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
