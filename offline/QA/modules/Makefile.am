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
	-lfun4all \
	-lphg4hit \
	-lg4detectors_io \
	-lcalo_io -lcalo_util \
	-lg4jets_io \
	-ltrackbase_historic_io \
	-lg4eval

pkginclude_HEADERS = \
	QAG4SimulationCalorimeter.h \
	QAG4SimulationCalorimeterSum.h \
	QAG4SimulationIntt.h \
	QAG4SimulationJet.h \
	QAG4SimulationMvtx.h \
	QAG4SimulationTpc.h \
	QAG4SimulationTracking.h \
	QAG4SimulationUpsilon.h \
	QAHistManagerDef.h

if ! MAKEROOT6
	ROOT5DICTS = \
		QAG4SimulationCalorimeter_Dict.cc \
		QAG4SimulationCalorimeterSum_Dict.cc \
		QAG4SimulationIntt_Dict.cc \
		QAG4SimulationJet_Dict.cc \
		QAG4SimulationMvtx_Dict.cc \
		QAG4SimulationTpc_Dict.cc \
		QAG4SimulationTracking_Dict.cc \
		QAG4SimulationUpsilon_Dict.cc \
		QAHistManagerDef_Dict.cc
endif

libqa_modules_la_SOURCES = \
	$(ROOT5DICTS) \
	QAG4SimulationCalorimeter.cc \
	QAG4SimulationCalorimeterSum.cc \
	QAG4SimulationIntt.cc \
	QAG4SimulationJet.cc \
	QAG4SimulationMvtx.cc \
	QAG4SimulationTpc.cc \
	QAG4SimulationTracking.cc \
	QAG4SimulationUpsilon.cc \
	QAHistManagerDef.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

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
	rm -f *Dict* testexternals.cc
