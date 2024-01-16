AUTOMAKE_OPTIONS = foreign

CINTINCLUDES = -pthread -I$(includedir) -I$(OFFLINE_MAIN)/lib -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include -I$(ROOTSYS)/lib/include -I$(OPT_SPHENIX)/include

INCLUDES = $(CINTINCLUDES)

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -lcalo_io \
  -lCLHEP 

pkginclude_HEADERS = \
  caloTowerEmbed.h

lib_LTLIBRARIES = \
  libCaloEmbedding.la

libCaloEmbedding_la_SOURCES = \
  caloTowerEmbed.cc \
  $(ROOTDICTS)

libCaloEmbedding_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libCaloembedding.la

CALOEMBEDDINGOBJS =  caloTowerEmbed.o derp.o

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f *Dict*  $(BUILT_SOURCES) *.pcm
