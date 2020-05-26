AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

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

if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
    ParticleFlowElement_Dict_rdict.pcm \
    ParticleFlowElementv1_Dict_rdict.pcm \
    ParticleFlowElementContainer_Dict_rdict.pcm 
endif

lib_LTLIBRARIES = \
  libparticleflow.la

libparticleflow_la_SOURCES = \
  $(ROOTDICTS) \
  ParticleFlowReco.cc \
  ParticleFlowElement.cc \
  ParticleFlowElementv1.cc \
  ParticleFlowElementContainer.cc \
  ParticleFlowJetInput.cc

libparticleflow_la_LIBADD = \
  -lg4dst \
  -lg4eval \
  -lcalo_io \
  -lphool \
  -lfun4all \
  -lSubsysReco

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

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
