AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) \
  -isystem $(OFFLINE_MAIN)/include/eigen3 \
  -isystem $(OFFLINE_MAIN)/include \
  -isystem $(G4_MAIN)/include \
  -isystem $(ROOTSYS)/include

RINCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I$(OFFLINE_MAIN)/include \
  -I$(G4_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

AM_CXXFLAGS =  -msse2

lib_LTLIBRARIES = libg4eval.la

libg4eval_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs` \
  -lfun4all \
  -lg4hough \
  -lHelixHough \
  -lFitNewton \
  -lcemc

pkginclude_HEADERS = \
  PHG4Evaluator.h \
  PHG4CalEvaluator.h

#pkginclude_HEADERS = $(include_HEADERS)

libg4eval_la_SOURCES = \
  EvalDict.C \
  PHG4Evaluator.C \
  PHG4CalEvaluator.C

EvalDict.C : \
  PHG4Evaluator.h \
  PHG4CalEvaluator.h \
  EvalLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(RINCLUDES) $^

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libg4eval.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
