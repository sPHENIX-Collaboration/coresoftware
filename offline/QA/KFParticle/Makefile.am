AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  QAKFParticle.h \
  QAKFParticleTrackPtAsymmetry.h

lib_LTLIBRARIES = \
  libkfparticleqa.la

libkfparticleqa_la_SOURCES = \
  QAKFParticle.cc \ 
  QAKFParticleTrackPtAsymmetry.cc

libkfparticleqa_io_la_LIBADD = \
  -lKFParticle \
  -lSubsysReco \
  -lqautils \
  -lphool


BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals 

testexternals_SOURCES = testexternals.cc

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD = libkfparticleqa.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)