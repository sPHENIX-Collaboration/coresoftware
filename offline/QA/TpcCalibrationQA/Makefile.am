AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = README.md

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

pkginclude_HEADERS = \
  TpcCalibrationQA.h \
  TpcConditionsReference.h \
  TpcPHGarfieldReference.h

noinst_HEADERS = \
  TpcQAInternal.h \
  TpcDistortionQA.h \
  TpcDriftQA.h \
  TpcConditionsUtils.h \
  TpcRunAverageConditions.h

lib_LTLIBRARIES = libtpccalibrationqa.la

libtpccalibrationqa_la_SOURCES = \
  TpcCalibrationQA.cc \
  TpcCalibrationQABatch.cc \
  TpcDistortionQA.cc \
  TpcDriftQA.cc \
  TpcQAInput.cc \
  TpcQAOutput.cc \
  TpcConditionsReference.cc \
  TpcPHGarfieldReference.cc \
  TpcRunAverageConditions.cc

libtpccalibrationqa_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lffamodules \
  -lcdbobjects \
  -lg4detectors_io \
  -lPHGarfield \
  -lTpcConditions \
  -lGraf3d

BUILT_SOURCES = testexternals.cc
noinst_PROGRAMS = testexternals
testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libtpccalibrationqa.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
