AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libQAKFParticleTracktoCalo.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

pkginclude_HEADERS = \
  QAKFParticleTracktoCalo.h

libQAKFParticleTracktoCalo_la_SOURCES = \
  QAKFParticleTracktoCalo.cc 

libQAKFParticleTracktoCalo_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lKFParticle \
  -lcalo_io \
  -lCLHEP \
  -lfun4all

clean-local:
	rm -f $(BUILT_SOURCES)