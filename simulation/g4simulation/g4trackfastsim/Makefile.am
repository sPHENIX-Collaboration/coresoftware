AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libg4trackfastsim.la 


AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(G4_MAIN)/include \
  -I$(OFFLINE_MAIN)/include  \
  -I$(OFFLINE_MAIN)/include/eigen3  \
  -I$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(ROOTSYS)/lib

libg4trackfastsim_la_LIBADD = \
  -lfun4all \
  -lphool \
  -lg4testbench \
  -lphg4hit \
  -ltrackbase_historic_io \
  -lcalo_io \
  -lgenfit2 \
  -lgenfit2exp \
  -lPHGenFit 

pkginclude_HEADERS = \
  PHG4TrackFastSim.h \
  PHG4TrackFastSimEval.h


if MAKEROOT6
  pcmdir = $(libdir)
else
  ROOT5_DICTS = \
    PHG4TrackFastSim_Dict.cc \
    PHG4TrackFastSimEval_Dict.cc
endif

libg4trackfastsim_la_SOURCES = \
  $(ROOT5_DICTS) \
  PHG4TrackFastSim.cc \
  PHG4TrackFastSimEval.cc



# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_tracker

testexternals_tracker_SOURCES = testexternals.cc
testexternals_tracker_LDADD = libg4trackfastsim.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
