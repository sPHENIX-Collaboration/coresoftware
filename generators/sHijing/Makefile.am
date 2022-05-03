AUTOMAKE_OPTIONS = foreign

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs-without-gui` \
  `fastjet-config --libs` \
  -lHepMC \
  -lhijing \
  -lhijing_dummy \
  -lgfortran \
  -lpthread

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include

AM_FFLAGS="-DFVOIDP=INTEGER*8"

noinst_HEADERS = \
  HijCrdn.h \
  HijJet1.h \
  HijJet2.h \
  HijJet4.h \
  HiMain1.h \
  HiMain2.h \
  HiParnt.h \
  HiStrng.h \
  RanSeed.h

bin_PROGRAMS = sHijing

noinst_PROGRAMS = xml_test

xml_test_SOURCES = xml_test.cc

sHijing_SOURCES = \
  Hijing.cxx \
  getaddr.cxx \
  ghijcrdn.F \
  ghijjet1.F \
  ghijjet2.F \
  ghijjet4.F \
  ghimain1.F \
  ghimain2.F \
  ghiparnt.F \
  ghistrng.F \
  granseed.F 


clean-local:
	rm -f $(BUILT_SOURCES)
