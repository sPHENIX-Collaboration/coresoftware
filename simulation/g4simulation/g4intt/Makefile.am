AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libg4intt_io.la \
    libg4intt.la

AM_CPPFLAGS = \
    -I$(includedir) \
    -isystem$(OFFLINE_MAIN)/include \
    -isystem$(ROOTSYS)/include

# set in configure.in to check if gcc version >= 4.8
#if GCC_GE_48
#  AM_CXXFLAGS = -std=c++11
#endif

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(ROOTSYS)/lib \
    -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  InttDeadMap.h \
  InttDeadMapv1.h \
  PHG4InttDeadMapLoader.h \
  PHG4InttDefs.h \
  PHG4InttDigitizer.h  \
  PHG4InttHitReco.h \
  PHG4InttSubsystem.h

ROOTDICTS = \
  InttDeadMap_Dict.cc \
  InttDeadMapv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  InttDeadMap_Dict_rdict.pcm \
  InttDeadMapv1_Dict_rdict.pcm

# sources for intt library
libg4intt_la_SOURCES = \
  PHG4InttDeadMapLoader.cc \
  PHG4InttDetector.cc \
  PHG4InttDigitizer.cc \
  PHG4InttDisplayAction.cc \
  PHG4InttFPHXParameterisation.cc \
  PHG4InttHitReco.cc \
  PHG4InttSteppingAction.cc \
  PHG4InttSubsystem.cc

libg4intt_la_LIBADD = \
  libg4intt_io.la \
  -lg4detectors \
  -lg4tracking_io \
  -lintt \
  -lintt_io

# sources for io library
libg4intt_io_la_SOURCES = \
  $(ROOTDICTS) \
  InttDeadMap.cc \
  InttDeadMapv1.cc

libg4intt_io_la_LIBADD = \
  -lphool \
  -lg4detectors_io 


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

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
