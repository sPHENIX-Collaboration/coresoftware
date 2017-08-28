AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libphfield_io.la \
   libphfield.la

AM_CXXFLAGS = \
 -Werror -msse2

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libphfield_io_la_LIBADD = \
  -lphool

libphfield_la_LIBADD = \
  -lfun4all \
  libphfield_io.la

pkginclude_HEADERS = \
  PHFieldConfig.h \
  PHFieldConfig_v1.h \
  PHFieldConfig_v2.h \
  PHFieldUtility.h \
  PHField.h

libphfield_io_la_SOURCES = \
  PHFieldConfig_Dict.C \
  PHFieldConfig.cc \
  PHFieldConfig_v1_Dict.C \
  PHFieldConfig_v1.cc \
  PHFieldConfig_v2_Dict.C \
  PHFieldConfig_v2.cc 

libphfield_la_SOURCES = \
	PHField.cc \
	PHFieldUniform.cc \
	PHField2D.cc \
	PHField3DCylindrical.cc \
	PHField3DCartesian.cc \
  PHFieldUtility_Dict.C \
	PHFieldUtility.cc 

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

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
	rm -f *Dict* $(BUILT_SOURCES)
