AUTOMAKE_OPTIONS = foreign


AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

# ThePEG headers are not clean under -Wextra -Wshadow -pedantic with gcc 14
#AM_CXXFLAGS = -Wno-shadow -Wno-unused-parameter -Wno-deprecated-declarations

lib_LTLIBRARIES = libPHHerwig7.la

pkginclude_HEADERS = \
  PHHerwig7.h \
  PHHerwig7Core.h \
  PHHerwig7GenTrigger.h \
  PHHerwig7JetTrigger.h

libPHHerwig7_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs` \
  `fastjet-config --libs`

libPHHerwig7_la_LIBADD = \
  -lSubsysReco \
  -lphhepmc \
  -lphool \
  -lHepMC \
  -lHerwigAPI \
  -lThePEG

libPHHerwig7_la_SOURCES = \
  PHHerwig7.cc \
  PHHerwig7GenTrigger.cc \
  PHHerwig7JetTrigger.cc

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHHerwig7.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
