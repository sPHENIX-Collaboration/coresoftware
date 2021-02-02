AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  FermimotionAfterburner.h

lib_LTLIBRARIES = \
  libFermimotionAfterburner.la

libFermimotionAfterburner_la_SOURCES = \
  FermimotionAfterburner.cc \
  FermiMotion.cc

libFermimotionAfterburner_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lphhepmc

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libFermimotionAfterburner.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
