AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = librhouecalib.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

librhouecalib_la_LIBADD = \
  -lcalo_io \
  -lcalotrigger_io \
  -lcdbobjects \
  -lglobalvertex_io \
  -ljetbackground_io \
  -lmbd \
  -lmbd_io \
  -lphool \
  -lSubsysReco

librhouecalib_la_SOURCES = \
  RhoUECalibBase.cc \
  MbdChargeHistMaker.cc \
  TowerShapeMapMaker.cc \
  RhoUECalibMaker.cc \
  RhoUECalibPayload.cc \
  RhoUECalibEventVeto.cc \
  MbdTTGainOverride.cc

pkginclude_HEADERS = \
  RhoUECalibBase.h \
  MbdChargeHistMaker.h \
  TowerShapeMapMaker.h \
  RhoUECalibMaker.h \
  RhoUECalibPayload.h \
  RhoUECalibEventVeto.h \
  MbdTTGainOverride.h

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  librhouecalib.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
