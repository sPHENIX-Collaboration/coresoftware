AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libtrack_reco.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include 
  
libtrack_reco_la_SOURCES = \
  PHInitVertexing.C \
  PHInitVertexing_Dict.C \
  PHTruthVertexing.C \
  PHTruthVertexing_Dict.C

libtrack_reco_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libtrack_reco_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lgsl \
  -lgslcblas \
  -lg4hough_io \
  -lphparameter
  
# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

##############################################
# please add new classes in alphabetical order

# pkginclude_HEADERS = 

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_track_reco


testexternals_track_reco_SOURCES = testexternals.C
testexternals_track_reco_LDADD = libtrack_reco.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
