AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libbcolumicount_io.la \
  libbcolumicount.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

libbcolumicount_io_la_LIBADD = \
  -lphool

libbcolumicount_la_LIBADD = \
  libbcolumicount_io.la \
  -lffaobjects \
  -lffarawobjects \
  -lSubsysReco

ROOTDICTS = \
  BcoInfo_Dict.cc \
  BcoInfov1_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

pkginclude_HEADERS = \
  BcoInfo.h \
  BcoInfov1.h \
  BcoLumiCheck.h \
  BcoLumiReco.h


libbcolumicount_io_la_SOURCES = \
  $(ROOTDICTS) \
  BcoInfo.cc \
  BcoInfov1.cc

libbcolumicount_la_SOURCES = \
  BcoLumiCheck.cc \
  BcoLumiReco.cc

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals \
  testexternals_io

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libbcolumicount.la

testexternals_io_SOURCES = \
  testexternals.cc

testexternals_io_LDADD = \
  libbcolumicount_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f $(BUILT_SOURCES)
