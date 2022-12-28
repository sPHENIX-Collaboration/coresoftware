##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libbbc.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I${G4_MAIN}/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  BbcOut.h \
  BbcOutV1.h \
  BbcNorthSouth.h \
  BbcNorthSouthV1.h \
  BbcPmtContainer.h \
  BbcPmtContainerV1.h \
  BbcPmtHit.h \
  BbcPmtHitV1.h \
  BbcReturnCodes.h

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  Bbc_Dict_rdict.pcm

libbbc_la_SOURCES = \
  Bbc_Dict.cc \
  BbcOut.cc \
  BbcOutV1.cc \
  BbcNorthSouth.cc \
  BbcNorthSouthV1.cc \
  BbcPmtHit.cc \
  BbcPmtHitV1.cc \
  BbcPmtContainer.cc \
  BbcPmtContainerV1.cc

libbbc_la_LIBADD = \
  -lphool

# Rule for generating table CINT dictionaries.
Bbc_Dict.cc: BbcOut.h  BbcOutV1.h BbcNorthSouth.h BbcNorthSouthV1.h BbcPmtContainer.h  BbcPmtContainerV1.h BbcPmtHit.h BbcPmtHitV1.h BbcLinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
Bbc_Dict_rdict.pcm: Bbc_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_bbc

testexternals_bbc_SOURCES = testexternals.cc
testexternals_bbc_LDADD = libbbc.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm

