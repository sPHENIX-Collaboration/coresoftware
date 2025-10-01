##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

lib_LTLIBRARIES = \
  libtpcchannelmasking.la

pkginclude_HEADERS = \
  TpcMaskedChannelMap.h

# sources for tpc library
libtpcchannelmasking_la_SOURCES = \
  TpcMaskedChannelMap.cc

libtpcchannelmasking_la_LIBADD = \
  -lffarawobjects \
  -lcdbobjects \
  -ltrack \
  -ltpc \
  -lfun4all \
  -lphool \
  -lpthread \
  -lg4detectors_io

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libtpcchannelmasking.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f $(BUILT_SOURCES)
