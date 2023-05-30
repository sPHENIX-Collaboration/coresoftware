##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libbbc_io.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I${G4_MAIN}/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  BbcDefs.h \
  BbcOut.h \
  BbcOutV1.h \
  BbcNorthSouth.h \
  BbcNorthSouthV1.h \
  BbcPmtContainer.h \
  BbcPmtContainerV1.h \
  BbcPmtHit.h \
  BbcPmtHitV1.h \
  BbcReconstruction.h \
  BbcReturnCodes.h \
  BbcVertex.h \
  BbcVertexv1.h \
  BbcVertexv2.h \
  BbcVertexMap.h \
  BbcVertexMapv1.h

ROOTDICTS = \
  BbcOut_Dict.cc \
  BbcOutV1_Dict.cc \
  BbcNorthSouth_Dict.cc \
  BbcNorthSouthV1_Dict.cc \
  BbcPmtHit_Dict.cc \
  BbcPmtHitV1_Dict.cc \
  BbcPmtContainer_Dict.cc \
  BbcPmtContainerV1_Dict.cc \
  BbcVertex_Dict.cc \
  BbcVertexv1_Dict.cc \
  BbcVertexv2_Dict.cc \
  BbcVertexMap_Dict.cc \
  BbcVertexMapv1_Dict.cc

pcmdir = $(libdir)

nobase_dist_pcm_DATA = \
  BbcOut_Dict_rdict.pcm \
  BbcOutV1_Dict_rdict.pcm \
  BbcNorthSouth_Dict_rdict.pcm \
  BbcNorthSouthV1_Dict_rdict.pcm \
  BbcPmtHit_Dict_rdict.pcm \
  BbcPmtHitV1_Dict_rdict.pcm \
  BbcPmtContainer_Dict_rdict.pcm \
  BbcPmtContainerV1_Dict_rdict.pcm \
  BbcVertex_Dict_rdict.pcm \
  BbcVertexv1_Dict_rdict.pcm \
  BbcVertexv2_Dict_rdict.pcm \
  BbcVertexMap_Dict_rdict.pcm \
  BbcVertexMapv1_Dict_rdict.pcm


libbbc_io_la_SOURCES = \
  $(ROOTDICTS) \
  BbcOut.cc \
  BbcOutV1.cc \
  BbcNorthSouth.cc \
  BbcNorthSouthV1.cc \
  BbcPmtHit.cc \
  BbcPmtHitV1.cc \
  BbcPmtContainer.cc \
  BbcPmtContainerV1.cc \
  BbcReconstruction.cc \
  BbcVertexv1.cc \
  BbcVertexv2.cc \
  BbcVertexMap.cc \
  BbcVertexMapv1.cc

libbbc_io_la_LIBADD = \
  -lphool \
  -lSubsysReco

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libbbc_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm

