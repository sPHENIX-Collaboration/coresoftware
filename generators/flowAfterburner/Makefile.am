AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libflowafterburner.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lCLHEP \
  -lgsl \
  -lgslcblas \
  -lHepMC

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include

pkginclude_HEADERS = \
  flowAfterburner.h

libflowafterburner_la_SOURCES = \
  flowAfterburner.cc


bin_PROGRAMS = flowAfterburner 

noinst_PROGRAMS = test

flowAfterburner_SOURCES = \
  main.cc

flowAfterburner_LDADD = \
  libflowafterburner.la

test_LDADD = \
  libflowafterburner.la

test_SOURCES = \
  test.cc

clean-local:
	rm -f $(BUILT_SOURCES)

