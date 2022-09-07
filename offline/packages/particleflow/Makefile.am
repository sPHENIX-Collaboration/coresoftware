AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  ParticleFlowReco.h \
  ParticleFlowElement.h \
  ParticleFlowElementv1.h \
  ParticleFlowElementContainer.h \
  ParticleFlowJetInput.h

ROOTDICTS = \
  ParticleFlowElement_Dict.cc \
  ParticleFlowElementv1_Dict.cc \
  ParticleFlowElementContainer_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  ParticleFlowElement_Dict_rdict.pcm \
  ParticleFlowElementv1_Dict_rdict.pcm \
  ParticleFlowElementContainer_Dict_rdict.pcm 

lib_LTLIBRARIES = \
  libparticleflow_io.la \
  libparticleflow.la

libparticleflow_io_la_SOURCES = \
  $(ROOTDICTS) \
  ParticleFlowElement.cc \
  ParticleFlowElementv1.cc \
  ParticleFlowElementContainer.cc

libparticleflow_la_SOURCES = \
  ParticleFlowReco.cc \
  ParticleFlowJetInput.cc

libparticleflow_io_la_LIBADD = \
  -lphool

libparticleflow_la_LIBADD = \
  libparticleflow_io.la \
  -lcalo_io \
  -lgsl \
  -lgslcblas \
  -lphg4hit \
  -lg4jets_io \
  -ltrackbase_historic_io \
  -lCLHEP \
  -lg4vertex_io \
  -lSubsysReco

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_io \
  testexternals

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD   = libparticleflow_io.la

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libparticleflow.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
