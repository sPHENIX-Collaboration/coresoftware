AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `fastjet-config --libs`

pkginclude_HEADERS = \
  HepMCJetTrigger.h

lib_LTLIBRARIES = \
  libHepMCJetTrigger.la

libHepMCJetTrigger_la_SOURCES = \
  HepMCJetTrigger.cc

libHepMCJetTrigger_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lphhepmc 

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libHepMCJetTrigger.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
