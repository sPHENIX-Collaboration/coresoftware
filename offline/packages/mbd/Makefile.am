##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
if USE_ONLINE
lib_LTLIBRARIES = \
  libmbd_io.la
else
lib_LTLIBRARIES = \
  libmbd_io.la \
  libmbd.la
endif

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I${G4_MAIN}/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(OFFLINE_MAIN)/lib

if USE_ONLINE
pkginclude_HEADERS = \
  BbcGeom.h \
  BbcGeomV1.h \
  MbdGeom.h \
  MbdGeomV1.h

else
pkginclude_HEADERS = \
  BbcDefs.h \
  MbdDefs.h \
  BbcOut.h \
  BbcOutV1.h \
  MbdOut.h \
  MbdOutV1.h \
  MbdOutV2.h \
  BbcPmtContainer.h \
  BbcPmtContainerV1.h \
  BbcPmtHit.h \
  BbcPmtHitV1.h \
  MbdPmtContainer.h \
  MbdPmtContainerV1.h \
  MbdPmtSimContainerV1.h \
  MbdPmtHit.h \
  MbdPmtHitV1.h \
  MbdPmtSimHitV1.h \
  BbcPmtInfoContainerV1.h \
  BbcPmtInfoV1.h \
  MbdSig.h \
  MbdEvent.h \
  MbdCalib.h \
  MbdReco.h \
  BbcGeom.h \
  BbcGeomV1.h \
  MbdGeom.h \
  MbdGeomV1.h \
  BbcReturnCodes.h \
  MbdReturnCodes.h \
  BbcVertex.h \
  BbcVertexv1.h \
  BbcVertexv2.h \
  BbcVertexMap.h \
  BbcVertexMapv1.h
endif

if USE_ONLINE
ROOTDICTS = \
  BbcGeom_Dict.cc \
  BbcGeomV1_Dict.cc \
  MbdGeom_Dict.cc \
  MbdGeomV1_Dict.cc

else
ROOTDICTS = \
  BbcOut_Dict.cc \
  BbcOutV1_Dict.cc \
  BbcGeom_Dict.cc \
  BbcGeomV1_Dict.cc \
  BbcPmtHit_Dict.cc \
  BbcPmtHitV1_Dict.cc \
  BbcPmtContainer_Dict.cc \
  BbcPmtContainerV1_Dict.cc \
  BbcPmtInfoV1_Dict.cc \
  BbcPmtInfoContainerV1_Dict.cc \
  BbcVertex_Dict.cc \
  BbcVertexv1_Dict.cc \
  BbcVertexv2_Dict.cc \
  BbcVertexMap_Dict.cc \
  BbcVertexMapv1_Dict.cc \
  MbdOut_Dict.cc \
  MbdOutV1_Dict.cc \
  MbdOutV2_Dict.cc \
  MbdGeom_Dict.cc \
  MbdGeomV1_Dict.cc \
  MbdPmtHit_Dict.cc \
  MbdPmtHitV1_Dict.cc \
  MbdPmtSimHitV1_Dict.cc \
  MbdPmtContainer_Dict.cc \
  MbdPmtContainerV1_Dict.cc \
  MbdPmtSimContainerV1_Dict.cc
endif

pcmdir = $(libdir)

if USE_ONLINE
nobase_dist_pcm_DATA = \
  BbcGeom_Dict_rdict.pcm \
  BbcGeomV1_Dict_rdict.pcm \
  MbdGeom_Dict_rdict.pcm \
  MbdGeomV1_Dict_rdict.pcm

else
nobase_dist_pcm_DATA = \
  BbcOut_Dict_rdict.pcm \
  BbcOutV1_Dict_rdict.pcm \
  BbcGeom_Dict_rdict.pcm \
  BbcGeomV1_Dict_rdict.pcm \
  BbcPmtHit_Dict_rdict.pcm \
  BbcPmtHitV1_Dict_rdict.pcm \
  BbcPmtContainer_Dict_rdict.pcm \
  BbcPmtContainerV1_Dict_rdict.pcm \
  BbcPmtInfoV1_Dict_rdict.pcm \
  BbcPmtInfoContainerV1_Dict_rdict.pcm \
  BbcVertex_Dict_rdict.pcm \
  BbcVertexv1_Dict_rdict.pcm \
  BbcVertexv2_Dict_rdict.pcm \
  BbcVertexMap_Dict_rdict.pcm \
  BbcVertexMapv1_Dict_rdict.pcm \
  MbdOut_Dict_rdict.pcm \
  MbdOutV1_Dict_rdict.pcm \
  MbdOutV2_Dict_rdict.pcm \
  MbdGeom_Dict_rdict.pcm \
  MbdGeomV1_Dict_rdict.pcm \
  MbdPmtHit_Dict_rdict.pcm \
  MbdPmtHitV1_Dict_rdict.pcm \
  MbdPmtSimHitV1_Dict_rdict.pcm \
  MbdPmtContainer_Dict_rdict.pcm \
  MbdPmtContainerV1_Dict_rdict.pcm \
  MbdPmtSimContainerV1_Dict_rdict.pcm
endif

if USE_ONLINE
libmbd_io_la_SOURCES = \
  $(ROOTDICTS) \
  BbcGeomV1.cc \
  MbdGeomV1.cc

else
libmbd_io_la_SOURCES = \
  $(ROOTDICTS) \
  BbcOut.cc \
  BbcOutV1.cc \
  BbcPmtHit.cc \
  BbcPmtHitV1.cc \
  BbcPmtContainer.cc \
  BbcPmtContainerV1.cc \
  BbcPmtInfoV1.cc \
  BbcPmtInfoContainerV1.cc \
  BbcGeomV1.cc \
  MbdOut.cc \
  MbdOutV1.cc \
  MbdOutV2.cc \
  MbdPmtHit.cc \
  MbdPmtHitV1.cc \
  MbdPmtSimHitV1.cc \
  MbdPmtContainer.cc \
  MbdPmtContainerV1.cc \
  MbdPmtSimContainerV1.cc \
  MbdGeomV1.cc \
  BbcVertexv1.cc \
  BbcVertexv2.cc \
  BbcVertexMap.cc \
  BbcVertexMapv1.cc

libmbd_la_SOURCES = \
  MbdSig.cc \
  MbdEvent.cc \
  MbdCalib.cc \
  MbdReco.cc

endif

if USE_ONLINE
libmbd_io_la_LIBADD = \
  -lphool \
  -lcdbobjects

else
libmbd_io_la_LIBADD = \
  -lcalo_io \
  -lphool

libmbd_la_LIBADD = \
  libmbd_io.la \
  -lphool \
  -lffamodules \
  -lcdbobjects \
  -lSubsysReco \
  -lglobalvertex_io

endif

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

if USE_ONLINE
noinst_PROGRAMS = \
  testexternals_io

else

noinst_PROGRAMS = \
  testexternals \
  testexternals_io

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libmbd.la

endif

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD = libmbd_io.la

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

