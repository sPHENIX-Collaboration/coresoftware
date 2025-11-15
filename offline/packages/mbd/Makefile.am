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

if USE_ONLINE
AM_CPPFLAGS = \
  -DONLINE \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I${G4_MAIN}/include \
  -isystem$(ROOTSYS)/include
else
AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I${G4_MAIN}/include \
  -isystem$(ROOTSYS)/include
endif

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(OFFLINE_MAIN)/lib

if USE_ONLINE
pkginclude_HEADERS = \
  MbdDefs.h \
  MbdEvent.h \
  MbdGeom.h \
  MbdGeomV1.h \
  MbdGeomV2.h \
  MbdOut.h \
  MbdOutV1.h \
  MbdOutV2.h \
  MbdPmtContainer.h \
  MbdPmtContainerV1.h \
  MbdPmtHit.h \
  MbdPmtHitV1.h \
  MbdRawContainer.h \
  MbdRawContainerV1.h \
  MbdRawHit.h \
  MbdRawHitV1.h \
  MbdReturnCodes.h \
  MbdRunningStats.h \
  MbdCalib.h \
  MbdSig.h

else
pkginclude_HEADERS = \
  MbdDefs.h \
  MbdOut.h \
  MbdOutV1.h \
  MbdOutV2.h \
  MbdPmtContainer.h \
  MbdPmtContainerV1.h \
  MbdPmtSimContainerV1.h \
  MbdPmtHit.h \
  MbdPmtHitV1.h \
  MbdPmtSimHitV1.h \
  MbdRawContainer.h \
  MbdRawContainerV1.h \
  MbdRawHit.h \
  MbdRawHitV1.h \
  MbdRunningStats.h \
  MbdSig.h \
  MbdEvent.h \
  MbdCalib.h \
  MbdReco.h \
  MbdGeom.h \
  MbdGeomV1.h \
  MbdGeomV2.h \
  MbdReturnCodes.h
endif

if USE_ONLINE
ROOTDICTS = \
  MbdGeom_Dict.cc \
  MbdGeomV1_Dict.cc \
  MbdGeomV2_Dict.cc \
  MbdOut_Dict.cc \
  MbdOutV1_Dict.cc \
  MbdOutV2_Dict.cc \
  MbdPmtContainer_Dict.cc \
  MbdPmtContainerV1_Dict.cc \
  MbdPmtHit_Dict.cc \
  MbdPmtHitV1_Dict.cc \
  MbdRawContainer_Dict.cc \
  MbdRawContainerV1_Dict.cc \
  MbdRawHit_Dict.cc \
  MbdRawHitV1_Dict.cc

else
ROOTDICTS = \
  MbdOut_Dict.cc \
  MbdOutV1_Dict.cc \
  MbdOutV2_Dict.cc \
  MbdGeom_Dict.cc \
  MbdGeomV1_Dict.cc \
  MbdGeomV2_Dict.cc \
  MbdPmtHit_Dict.cc \
  MbdPmtHitV1_Dict.cc \
  MbdPmtSimHitV1_Dict.cc \
  MbdPmtContainer_Dict.cc \
  MbdPmtContainerV1_Dict.cc \
  MbdPmtSimContainerV1_Dict.cc \
  MbdRawHit_Dict.cc \
  MbdRawHitV1_Dict.cc \
  MbdRawContainer_Dict.cc \
  MbdRawContainerV1_Dict.cc
endif

pcmdir = $(libdir)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

if USE_ONLINE
libmbd_io_la_SOURCES = \
  $(ROOTDICTS) \
  MbdEvent.cc \
  MbdGeomV1.cc \
  MbdGeomV2.cc \
  MbdOut.cc \
  MbdOutV1.cc \
  MbdOutV2.cc \
  MbdPmtHit.cc \
  MbdPmtHitV1.cc \
  MbdPmtContainer.cc \
  MbdPmtContainerV1.cc \
  MbdRawHit.cc \
  MbdRawHitV1.cc \
  MbdRawContainer.cc \
  MbdRawContainerV1.cc \
  MbdRunningStats.cc \
  MbdCalib.cc \
  MbdSig.cc

else
libmbd_io_la_SOURCES = \
  $(ROOTDICTS) \
  MbdGeomV1.cc \
  MbdGeomV2.cc \
  MbdOut.cc \
  MbdOutV1.cc \
  MbdOutV2.cc \
  MbdPmtHit.cc \
  MbdPmtHitV1.cc \
  MbdPmtSimHitV1.cc \
  MbdPmtContainer.cc \
  MbdPmtContainerV1.cc \
  MbdPmtSimContainerV1.cc \
  MbdRawHit.cc \
  MbdRawHitV1.cc \
  MbdRawContainer.cc \
  MbdRawContainerV1.cc \
  MbdRunningStats.cc \
  MbdSig.cc

libmbd_la_SOURCES = \
  MbdEvent.cc \
  MbdCalib.cc \
  MbdReco.cc \
  MbdRunningStats.cc \
  MbdSig.cc

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
  -lcdbobjects \
  -lffamodules \
  -lffaobjects \
  -lffarawobjects \
  -lglobalvertex_io \
  -lphool \
  -lphg4hit \
  -lSubsysReco

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

