AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libzdcinfo_io.la \
   libzdcinfo.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libzdcinfo_io_la_LIBADD = \
  -lphool

libzdcinfo_la_LIBADD = \
  libzdcinfo_io.la \
  -lfun4all \
  -lcalo_io \
  -lcdbobjects \
  -lffamodules \
  -lglobalvertex_io

pkginclude_HEADERS = \
  Zdcinfo.h \
  Zdcinfov1.h \
  Zdcinfov2.h \
  ZdcReco.h

ROOTDICTS = \
  Zdcinfo_Dict.cc \
  Zdcinfov1_Dict.cc \
  Zdcinfov2_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libzdcinfo_io_la_SOURCES = \
  $(ROOTDICTS) \
  Zdcinfo.cc \
  Zdcinfov1.cc \
  Zdcinfov2.cc


libzdcinfo_la_SOURCES = \
  ZdcReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_zdcinfo_io \
  testexternals_zdcinfo

BUILT_SOURCES = \
  testexternals.cc

testexternals_zdcinfo_io_SOURCES = testexternals.cc
testexternals_zdcinfo_io_LDADD = libzdcinfo_io.la

testexternals_zdcinfo_SOURCES = testexternals.cc
testexternals_zdcinfo_LDADD = libzdcinfo.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
