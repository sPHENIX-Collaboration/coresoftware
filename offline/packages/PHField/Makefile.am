AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(G4_MAIN)/include \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libphfield_io.la \
   libphfield.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libphfield_io_la_LIBADD = \
  -lphool

libphfield_la_LIBADD = \
  libphfield_io.la \
  -lfun4all

pkginclude_HEADERS = \
  PHField3DCartesian.h \
  PHFieldConfig.h \
  PHFieldConfigv1.h \
  PHFieldConfigv2.h \
  PHFieldUtility.h \
  PHField.h

ROOTDICTS = \
  PHFieldConfig_Dict.cc \
  PHFieldConfigv1_Dict.cc \
  PHFieldConfigv2_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  PHFieldConfig_Dict_rdict.pcm \
  PHFieldConfigv1_Dict_rdict.pcm \
  PHFieldConfigv2_Dict_rdict.pcm

libphfield_io_la_SOURCES = \
  $(ROOTDICTS) \
  PHFieldConfig.cc \
  PHFieldConfigv1.cc \
  PHFieldConfigv2.cc 

libphfield_la_SOURCES = \
  PHFieldUniform.cc \
  PHField2D.cc \
  PHField3DCylindrical.cc \
  PHField3DCartesian.cc \
  PHFieldUtility.cc 

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.C

noinst_PROGRAMS = \
  testexternals_phfield_io \
  testexternals_phfield


testexternals_phfield_io_SOURCES = testexternals.C
testexternals_phfield_io_LDADD = libphfield_io.la

testexternals_phfield_SOURCES = testexternals.C
testexternals_phfield_LDADD = libphfield.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
