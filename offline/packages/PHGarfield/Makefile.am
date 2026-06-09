AUTOMAKE_OPTIONS = foreign

ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs)

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  PHGarfield.h

lib_LTLIBRARIES = \
  libPHGarfield.la

libPHGarfield_la_LIBADD = \
  -lffamodules \
  -lffarawobjects \
  -lcdbobjects \
  -ltpc \
  -lNoRootEvent \
  -lphool \
  -lphfield \
  -lGarfield \
  -lSubsysReco

ROOTDICTS = \
  PHGarfield_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  PHGarfield_Dict_rdict.pcm

libPHGarfield_la_SOURCES = \
  $(ROOTDICTS) \
  PHGarfield.cc

bin_PROGRAMS = \
  GasModel \
  MergeGasFiles

MergeGasFiles_SOURCES = MergeGasFiles.cc

MergeGasFiles_LDFLAGS = \
  -L$(OFFLINE_MAIN)/lib64 \
  $(ROOT_LIBS)

MergeGasFiles_LDADD = \
  -lGarfield

GasModel_SOURCES = GasModel.cc

GasModel_LDFLAGS = \
  -L$(OFFLINE_MAIN)/lib64 \
  $(ROOT_LIBS)

GasModel_LDADD = \
  -lGarfield

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libPHGarfield.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
