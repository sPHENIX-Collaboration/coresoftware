AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libcalotemp.la

bin_PROGRAMS = testODBC

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include

pkginclude_HEADERS = \
  CaloTemp.h

libcalotemp_la_SOURCES = \
  CaloTemp.cc 

testODBC_SOURCES = testODBC.cc
testODBC_LDADD = libcalotemp.la

libcalotemp_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  -lcalo_io \
  -lodbc++ \
  -lfun4all \
  -lffamodules

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libcalotemp.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
