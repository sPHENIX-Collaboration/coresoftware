AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(top_srcdir)/helix_hough/Kalman \
  -I$(top_srcdir)/Seamstress \
  -I$(top_srcdir)/FitNewton \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OFFLINE_MAIN)/include/eigen3

AM_CXXFLAGS = -DNDEBUG 

lib_LTLIBRARIES = \
  libHelixHough.la

noinst_HEADERS = \
  vector_math_inline.h

pkginclude_HEADERS = \
  HelixHough.h \
  fastvec.h \
  HelixResolution.h \
  HelixRange.h \
  SimpleHit3D.h \
  SimpleTrack3D.h \
  VertexFinder.h \
  VertexFitFunc.h \
  Kalman/HelixKalman.h \
  Kalman/HelixKalmanState.h \
  Kalman/CylinderKalman.h \
  sPHENIX/sPHENIXSeedFinder.h \
  sPHENIX/sPHENIXTracker.h \
  sPHENIX/sPHENIXTrackerTPC.h

libHelixHough_la_SOURCES = \
  SimpleHit3D.cpp \
  sPHENIX/sPHENIXSeedFinder.cpp \
  sPHENIX/sPHENIXTracker.cpp \
  sPHENIX/sPHENIXTrackerTPC.cpp \
  HelixHough_findHelices.cpp \
  HelixHough_findPairs.cpp \
  HelixHough_init.cpp \
  HelixHough_phiRange_sse.cpp \
  HelixHough_vote_sse.cpp \
  HelixHough_vote_pairs_sse.cpp \
  HelixHough_allButKappaRange_sse.cpp \
  HelixHough_dzdlRange_sse.cpp \
  HelixHough_split.cpp \
  VertexFinder.cpp \
  VertexFitFunc.cpp \
  Kalman/HelixKalman.cpp \
  Kalman/HelixKalmanState.cpp \
  Kalman/CylinderKalman.cpp


libHelixHough_la_LIBADD = \
$(top_builddir)/Seamstress/libSeamstress.la \
$(top_builddir)/FitNewton/libFitNewton.la


################################################
# linking tests

noinst_PROGRAMS = \
  testexternals 

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libHelixHough.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)




