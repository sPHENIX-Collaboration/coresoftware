AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libg4intt_io.la \
    libg4intt.la

AM_CPPFLAGS = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(ROOTSYS)/include \
    -I$(G4_MAIN)/include

# set in configure.in to check if gcc version >= 4.8
#if GCC_GE_48
#  AM_CXXFLAGS = -std=c++11
#endif

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libg4intt_io_la_LIBADD = \
  -lphool \
  -lg4detectors_io 

libg4intt_la_LIBADD = \
  libg4intt_io.la \
  -lg4detectors \
  -lintt_io

pkginclude_HEADERS = \
  INTTDeadMap.h \
  INTTDeadMapv1.h \
  PHG4InttDigitizer.h  \
  PHG4INTTDefs.h \
  PHG4INTTHitReco.h \
  PHG4INTTSubsystem.h

ROOTDICTS = \
  INTTDeadMap_Dict.cc \
  INTTDeadMapv1_Dict.cc
if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
  INTTDeadMap_Dict_rdict.pcm \
  INTTDeadMapv1_Dict_rdict.pcm
else
  ROOT5_DICTS = \
    PHG4INTTDeadMapLoader_Dict.cc \
    PHG4INTTDefs_Dict.cc \
    PHG4InttDigitizer_Dict.cc \
    PHG4INTTHitReco_Dict.cc \
    PHG4INTTSubsystem_Dict.cc
endif

libg4intt_io_la_SOURCES = \
  $(ROOTDICTS) \
  INTTDeadMap.cc \
  INTTDeadMapv1.cc

libg4intt_la_SOURCES = \
  $(ROOT5_DICTS) \
  PHG4INTTDeadMapLoader.cc \
  PHG4INTTDetector.cc \
  PHG4INTTDisplayAction.cc \
  PHG4InttDigitizer.cc \
  PHG4INTTHitReco.cc \
  PHG4INTTFPHXParameterisation.cc \
  PHG4INTTSteppingAction.cc \
  PHG4INTTSubsystem.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_g4intt \
  testexternals_g4intt_io

testexternals_g4intt_SOURCES = testexternals.cc
testexternals_g4intt_LDADD = libg4intt.la

testexternals_g4intt_io_SOURCES = testexternals.cc
testexternals_g4intt_io_LDADD = libg4intt_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
