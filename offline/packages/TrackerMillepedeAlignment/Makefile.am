##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libtrackeralign.la

libtrackeralign_la_LIBADD = \
  -lFROG \
  -lphool \
  -lSubsysReco \
  -lg4detectors_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -ltpc_io

pkginclude_HEADERS = \
  MakeMilleFiles.h \
  Mille.h \
  HelicalFitter.h

pcmdir = $(libdir)

libtrackeralign_la_SOURCES = \
  MakeMilleFiles.cc \
  Mille.cc \
  HelicalFitter.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_trackeralign

testexternals_trackeralign_SOURCES = testexternals.cc
testexternals_trackeralign_LDADD = libtrackeralign.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f $(BUILT_SOURCES)
