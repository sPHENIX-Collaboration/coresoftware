AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libqa_modules.la

libqa_modules_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libqa_modules_la_LIBADD = \
  -lcalo_io \
  -lcalo_util \
  -lfun4all \
  -lg4detectors_io \
  -lg4eval \
  -lg4jets_io \
  -lphg4hit \
  -ltrackbase_historic_io

pkginclude_HEADERS = \
  QAG4SimulationCalorimeter.h \
  QAG4SimulationCalorimeterSum.h \
  QAG4SimulationIntt.h \
  QAG4SimulationJet.h \
  QAG4SimulationMvtx.h \
  QAG4SimulationTpc.h \
  QAG4SimulationTracking.h \
  QAG4SimulationVertex.h \
  QAG4SimulationUpsilon.h \
  QAHistManagerDef.h

libqa_modules_la_SOURCES = \
  QAG4SimulationCalorimeter.cc \
  QAG4SimulationCalorimeterSum.cc \
  QAG4SimulationIntt.cc \
  QAG4SimulationJet.cc \
  QAG4SimulationMvtx.cc \
  QAG4SimulationTpc.cc \
  QAG4SimulationTracking.cc \
  QAG4SimulationVertex.cc \
  QAG4SimulationUpsilon.cc \
  QAHistManagerDef.cc

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
	testexternals.cc

testexternals_SOURCES = \
	testexternals.cc

testexternals_LDADD = \
	libqa_modules.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
