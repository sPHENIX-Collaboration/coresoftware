AM_CXXFLAGS = -I$(top_srcdir)/helix_hough -I$(OFFLINE_MAIN)/include/eigen3 -I$(top_srcdir)/Seamstress -I$(top_srcdir)/FitNewton -I$(OFFLINE_MAIN)/include `root-config --cflags`

AM_LDFLAGS = ../libHelixHough.la ../../FitNewton/libFitNewton.la @ROOTLIBS@

test_with_vertex_SOURCES = test_with_vertex.cpp

bin_PROGRAMS = test_with_vertex
