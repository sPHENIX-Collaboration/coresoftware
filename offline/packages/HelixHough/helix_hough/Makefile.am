AM_CXXFLAGS = -msse2 -O2 -DNDEBUG -I$(OFFLINE_MAIN)/include/eigen3 -I$(top_srcdir)/Seamstress -I$(top_srcdir)/FitNewton -I$(OFFLINE_MAIN)/include
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
$(eigen3_LIBS) \
$(fitnewton_LIBS) \
$(seamstress_LIBS)

libHelixHough_la_LDFLAGS = -version-info 1:0:0

libHelixHough_la_CPPFLAGS = \
-I$(top_srcdir)/helix_hough/sPHENIX \
-I$(top_srcdir)/helix_hough/Kalman \
$(eigen3_CFLAGS) \
$(fitnewton_CFLAGS)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = helix_hough.pc





