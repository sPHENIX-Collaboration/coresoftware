AUTOMAKE_OPTIONS = foreign

INCLUDES = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
    libTPCbase.la

AM_CPPFLAGS = -DCGAL_DISABLE_ROUNDING_MATH_CHECK

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libTPCbase_la_LIBADD = \
  -L$(ROOTSYS)/lib \
  -lphool  \
  -lboost_filesystem \
  -lboost_system \
  -lXMLIO \
  -lCGAL

pkginclude_HEADERS = \
  TPCDataTypes.h \
  TPCConstants.h \
  TPCPadMap.h \
  vHit.h \
  TPCHit.h \
  TPCDigit.h \
  vCluster.h \
  TPCHitsContainer.h \
  TPCDigitsContainer.h

libTPCbase_la_SOURCES = \
  TPCPadMap.cxx \
  vHit.cxx \
  TPCHit.cxx \
  TPCDigit.cxx \
  vCluster.cxx \
  TPCHitsContainer.cxx \
  TPCDigitsContainer.cxx \
  TPCbase_Dict.cc

# Rule for generating table CINT dictionaries.
TPCbase_Dict.cc: TPCPadMap.h TPCbaseLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_TPCbase
testexternals_TPCbase_SOURCES = testexternals.cc
testexternals_TPCbase_LDADD = libTPCbase.la

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
