AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -DG4INTY_USE_XT -DG4VIS_USE_OPENGL -DG4UI_USE_TCSH -DG4UI_USE_XM -DG4VIS_USE_OPENGLXM -DG4VIS_USE_RAYTRACERX -DG4VIS_USE_OPENGLX

# set in configure.in to check if gcc version >= 4.8
#if GCC_GE_48
#  AM_CXXFLAGS += -std=c++11
#endif

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libphg4hit.la \
  libg4testbench.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(XERCESCROOT)/include \
  -isystem $(OFFLINE_MAIN)/include \
  -isystem $(ROOTSYS)/include \
  -isystem ${G4_MAIN}/include \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I$(OPT_SPHENIX)/include

# rootcint barfs on isystem
RINCLUDES = \
  -I$(includedir) \
  -I$(XERCESCROOT)/include \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I${G4_MAIN}/include \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I$(OPT_SPHENIX)/include


libphg4hit_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib


libphg4hit_la_LIBADD = \
  -lphool

libg4testbench_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs`

libg4testbench_la_LIBADD = \
  libphg4hit.la \
  -lboost_filesystem \
  -lEICPhysicsList \
  -leicsmear \
  -lg4decayer \
  -lgsl \
  -lgslcblas \
  -lHepMC \
  -lphfield \
  -lphgeom \
  -lphg4gdml \
  -lphhepmc \
  -lvararray

# I/O dictionaries have to exist for root5 and root6. For ROOT6 we need
# pcm files in addition. If someone can figure out how to make a list
# so this list of dictionaries is transformed into a list of pcm files
# following a simple naming convention, please change this accordingly and
# let me know
ROOTHITDICTS = \
  PHG4EventHeader_Dict.cc \
  PHG4EventHeaderv1_Dict.cc \
  PHG4Hit_Dict.cc \
  PHG4Hitv1_Dict.cc \
  PHG4HitEval_Dict.cc \
  PHG4HitContainer_Dict.cc \
  PHG4InEvent_Dict.cc \
  PHG4Particle_Dict.cc \
  PHG4Particlev1_Dict.cc \
  PHG4Particlev2_Dict.cc \
  PHG4Particlev3_Dict.cc \
  PHG4Shower_Dict.cc \
  PHG4Showerv1_Dict.cc \
  PHG4VtxPoint_Dict.cc \
  PHG4VtxPointv1_Dict.cc \
  PHG4TruthInfoContainer_Dict.cc

# for root6 we need pcm and dictionaries but only for
# i/o classes. For root5 we need only dictionaries but
# those for i/o and classes available on the cmd line
# MAKEROOT6 is set in the configure.ac, true for root6
if MAKEROOT6
# this is a tweak to install files in $(libdir), automake refuses
# to install other files in libdir, this construct gets around this
pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  PHG4EventHeader_Dict_rdict.pcm \
  PHG4EventHeaderv1_Dict_rdict.pcm \
  PHG4Hit_Dict_rdict.pcm \
  PHG4Hitv1_Dict_rdict.pcm \
  PHG4HitEval_Dict_rdict.pcm \
  PHG4HitContainer_Dict_rdict.pcm \
  PHG4InEvent_Dict_rdict.pcm \
  PHG4Particle_Dict_rdict.pcm \
  PHG4Particlev1_Dict_rdict.pcm \
  PHG4Particlev2_Dict_rdict.pcm \
  PHG4Particlev3_Dict_rdict.pcm \
  PHG4Shower_Dict_rdict.pcm \
  PHG4Showerv1_Dict_rdict.pcm \
  PHG4VtxPoint_Dict_rdict.pcm \
  PHG4VtxPointv1_Dict_rdict.pcm \
  PHG4TruthInfoContainer_Dict_rdict.pcm
else
  ROOT5HITDICTS = PHG4HitDefs_Dict.cc
  ROOT5TBDICTS = \
    HepMCNodeReader_Dict.cc \
    PHG4ConsistencyCheck_Dict.cc \
    PHG4HeadReco_Dict.cc \
    PHG4InEventCompress_Dict.cc \
    PHG4InEventReadBack_Dict.cc \
    PHG4InputFilter_Dict.cc \
    PHG4IonGun_Dict.cc \
    PHG4ParticleGenerator_Dict.cc \
    PHG4ParticleGeneratorBase_Dict.cc \
    PHG4ParticleGeneratorVectorMeson_Dict.cc \
    PHG4ParticleGeneratorD0_Dict.cc \
    PHG4ParticleGun_Dict.cc \
    PHG4Reco_Dict.cc \
    PHG4ScoringManager_Dict.cc \
    PHG4SimpleEventGenerator_Dict.cc \
    PHG4PileupGenerator_Dict.cc \
    PHG4Subsystem_Dict.cc \
    PHG4TruthSubsystem_Dict.cc \
    PHG4Utils_Dict.cc \
    ReadEICFiles_Dict.cc \
    PHG4HitReadBack_Dict.cc
endif

libphg4hit_la_SOURCES = \
  $(ROOTHITDICTS) \
  $(ROOT5HITDICTS) \
  PHG4EventHeaderv1.cc \
  PHG4Hit.cc \
  PHG4Hitv1.cc \
  PHG4HitContainer.cc \
  PHG4HitDefs.cc \
  PHG4HitEval.cc \
  PHG4InEvent.cc \
  PHG4Particle.cc \
  PHG4Particlev1.cc \
  PHG4Particlev2.cc \
  PHG4Particlev3.cc \
  PHG4Showerv1.cc \
  PHG4TruthInfoContainer.cc \
  PHG4VtxPoint.cc \
  PHG4VtxPointv1.cc

libg4testbench_la_SOURCES = \
  $(ROOT5TBDICTS) \
  G4TBMagneticFieldSetup.cc \
  G4TBFieldMessenger.cc \
  HepMCNodeReader.cc \
  PHG4ConsistencyCheck.cc \
  PHG4DisplayAction.cc \
  PHG4Detector.cc \
  PHG4HitReadBack.cc \
  PHG4EtaParameterization.cc \
  PHG4EtaPhiParameterization.cc \
  PHG4HeadReco.cc \
  PHG4InEventCompress.cc \
  PHG4InEventReadBack.cc \
  PHG4InputFilter.cc \
  PHG4IonGun.cc \
  PHG4MagneticField.cc \
  PHG4ParticleGeneratorBase.cc \
  PHG4ParticleGenerator.cc \
  PHG4ParticleGeneratorVectorMeson.cc \
  PHG4ParticleGeneratorD0.cc \
  PHG4ParticleGun.cc \
  PHG4PhenixDetector.cc \
  PHG4PhenixDisplayAction.cc \
  PHG4PhenixEventAction.cc \
  PHG4PhenixSteppingAction.cc \
  PHG4PhenixTrackingAction.cc \
  PHG4PileupGenerator.cc \
  PHG4PrimaryGeneratorAction.cc \
  PHG4Reco.cc \
  PHG4RegionInformation.cc \
  PHG4ScoringManager.cc \
  PHG4SimpleEventGenerator.cc \
  PHG4SteppingAction.cc \
  PHG4TrackUserInfoV1.cc \
  PHG4TruthEventAction.cc \
  PHG4TruthSubsystem.cc \
  PHG4TruthTrackingAction.cc \
  PHG4UIsession.cc \
  PHG4Utils.cc \
  ReadEICFiles.cc


##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  HepMCNodeReader.h \
  PHBBox.h \
  PHG4ColorDefs.h \
  PHG4Detector.h \
  PHG4DisplayAction.h \
  PHG4EventAction.h \
  PHG4EventHeader.h \
  PHG4HitDefs.h \
  PHG4Hit.h \
  PHG4Hitv1.h \
  PHG4HitEval.h \
  PHG4HitContainer.h \
  PHG4InEvent.h \
  PHG4IonGun.h \
  PHG4Particle.h \
  PHG4Particlev1.h \
  PHG4Particlev2.h \
  PHG4ParticleGenerator.h \
  PHG4ParticleGeneratorBase.h \
  PHG4ParticleGeneratorVectorMeson.h \
  PHG4ParticleGun.h \
  PHG4PhenixDetector.h \
  PHG4PileupGenerator.h \
  PHG4PrimaryGeneratorAction.h \
  PHG4Reco.h \
  PHG4RegionInformation.h \
  PHG4SimpleEventGenerator.h \
  PHG4Shower.h \
  PHG4Showerv1.h \
  PHG4SteppingAction.h \
  PHG4Subsystem.h \
  PHG4TrackingAction.h \
  PHG4TrackUserInfoV1.h \
  PHG4TruthInfoContainer.h \
  PHG4TruthSubsystem.h \
  PHG4Units.h \
  PHG4Utils.h \
  PHG4VtxPoint.h \
  PHG4VtxPointv1.h \
  ReadEICFiles.h

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4hits \
  testexternals_g4tb

BUILT_SOURCES = testexternals.cc

testexternals_g4hits_SOURCES = testexternals.cc
testexternals_g4hits_LDADD = libphg4hit.la

testexternals_g4tb_SOURCES = testexternals.cc
testexternals_g4tb_LDADD = libg4testbench.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(RINCLUDES) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm 
