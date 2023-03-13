AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = libfieldsim.la   

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -lgfortran \
  -lphool \
  -lSubsysReco

libfieldsim_la_SOURCES = \
  AnnularFieldSim.cc \
  AnalyticFieldModel.cc \
  ChargeMapReader.cc \
  Rossegger.cc \
  src.f \
  airy.f \
  d1mach.f 

pkginclude_HEADERS = \
  AnnularFieldSim.h \
  AnalyticFieldModel.h \
  ChargeMapReader.h \
  MultiArray.h \
  Rossegger.h

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libfieldsim.la

testexternals.cc:
	echo "this is a test.  $(ROOTSYS)/include"
	echo $(ROOTSYS)
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
