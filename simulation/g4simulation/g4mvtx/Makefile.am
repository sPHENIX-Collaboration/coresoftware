AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libg4mvtx_io.la \
  libg4mvtx.la 


AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(G4_MAIN)/include \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

libg4mvtx_io_la_LIBADD = \
  -lg4detectors_io 

libg4mvtx_la_LIBADD = \
  libg4mvtx_io.la \
  -lfun4all \
  -lphool \
  -lg4testbench \
  -lphg4hit \
  -lg4detectors \
  -ltrackbase_historic_io \
  -ltrack_io \
  -lmvtx

pkginclude_HEADERS = \
  PHG4MVTXHitReco.h \
  PHG4MVTXSubsystem.h \
  PHG4CylinderCell_MVTX.h \
  PHG4CylinderGeom_MVTX.h \
  PHG4MVTXDigitizer.h 

ROOTDICTS = \
  PHG4CylinderCell_MVTX_Dict.cc \
  PHG4CylinderGeom_MVTX_Dict.cc
if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
    PHG4CylinderCell_MVTX_Dict_rdict.pcm \
    PHG4CylinderGeom_MVTX_Dict_rdict.pcm
else
  ROOT5_DICTS = \
    PHG4MVTXCellReco_Dict.cc \
    PHG4MVTXSubsystem_Dict.cc \
    PHG4MVTXDigitizer_Dict.cc 
endif

libg4mvtx_la_SOURCES = \
  PHG4MVTXHitReco.cc \
  PHG4MVTXHitReco_Dict.cc \
  $(ROOT5_DICTS) \
  PHG4MVTXCellReco.cc \
  PHG4MVTXDetector.cc \
  PHG4MVTXSteppingAction.cc \
  PHG4MVTXSubsystem.cc \
  PHG4CylinderGeom_MVTX.cc \
  PHG4MVTXDigitizer.cc


libg4mvtx_io_la_SOURCES = \
  $(ROOTDICTS) \
  PHG4CylinderCell_MVTX.cc


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
testexternals_tracker_LDADD = libg4mvtx.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
