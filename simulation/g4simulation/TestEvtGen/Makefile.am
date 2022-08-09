AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libevtgendecayer.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -DEVTGEN_HEPMC3

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -I$(OFFLINE_MAIN)/include \
  `root-config --libs` 


libevtgendecayer_la_LIBADD = \
  -lfun4all \
  -lg4eval \
  -lEvtGen \
  -lEvtGenExternal \
  -lHepMC3

libevtgendecayer_la_SOURCES = \
  PHEvtGenDecayer.cxx 


pkginclude_HEADERS = \
  PHEvtGenDecayer.h 


################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libevtgendecayer.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
