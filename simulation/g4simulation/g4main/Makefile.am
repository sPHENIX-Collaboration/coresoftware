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
  -L$(OFFLINE_MAIN)/lib \
  -L$(XERCESCROOT)/lib \
  -L$(OPT_SPHENIX)/lib \
  `geant4-config --libs`


libphg4hit_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -leicsmear

libphg4hit_la_SOURCES = \
    PHG4AllHit_Dict.cc \
    PHG4EventHeader.cc \
    PHG4EventHeaderv1.cc \
    PHG4HitReadBack.cc \
    PHG4HitDefs.cc \
    PHG4Hit.cc \
    PHG4Hit_Dict.cc \
    PHG4Hitv1.cc \
    PHG4HitEval.cc \
    PHG4HitContainer.cc \
    PHG4InEvent.cc \
    PHG4Particle.cc \
    PHG4Particle_Dict.cc \
    PHG4Particlev1.cc \
    PHG4Particlev2.cc \
    PHG4Shower.cc \
    PHG4Showerv1.cc \
    PHG4TruthInfoContainer.cc \
    PHG4VtxPoint.cc \
    PHG4VtxPointv1.cc

libg4testbench_la_SOURCES = \
    G4TBMagneticFieldSetup.cc \
    G4TBFieldMessenger.cc \
    PHG4MagneticField.cc \
    HepMCNodeReader.cc \
    PHG4ConsistencyCheck.cc \
    PHG4Detector.cc \
    PHG4EtaParameterization.cc \
    PHG4EtaPhiParameterization.cc \
    PHG4HeadReco.cc \
    PHG4InEventCompress.cc \
    PHG4InEventReadBack.cc \
    PHG4InputFilter.cc \
    PHG4PileupGenerator.cc \
    PHG4SimpleEventGenerator.cc \
    PHG4ParticleGun.cc \
    PHG4ParticleGeneratorBase.cc \
    PHG4ParticleGenerator.cc \
    PHG4ParticleGeneratorVectorMeson.cc \
    PHG4ParticleGeneratorD0.cc \
    PHG4PhenixDetector.cc \
    PHG4PhenixEventAction.cc \
    PHG4PhenixSteppingAction.cc \
    PHG4PhenixTrackingAction.cc \
    PHG4_Dict.cc \
    PHG4PrimaryGeneratorAction.cc \
    PHG4Reco.cc \
    PHG4ScoringManager.cc \
    PHG4RegionInformation.cc \
    PHG4TrackUserInfoV1.cc \
    PHG4TruthEventAction.cc \
    PHG4TruthSteppingAction.cc \
    PHG4TruthSubsystem.cc \
    PHG4TruthTrackingAction.cc \
    PHG4SteppingAction.cc \
    PHG4UIsession.cc \
    PHG4Utils.cc \
    ReadEICFiles.cc

libg4testbench_la_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib \
    -lboost_filesystem \
    -lphfield \
    -lg4decayer \
    -lgsl -lgslcblas \
    -lvararray \
    -lHepMC \
    -lphhepmc \
    -lxerces-c \
    `geant4-config --libs` \
    -lphgeom \
    -lphg4gdml

libg4testbench_la_LIBADD = \
    libphg4hit.la

##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  HepMCNodeReader.h \
  PHBBox.h \
  PHG4ColorDefs.h \
  PHG4Detector.h \
  PHG4EventAction.h \
  PHG4EventHeader.h \
  PHG4HitDefs.h \
  PHG4Hit.h \
  PHG4Hitv1.h \
  PHG4HitEval.h \
  PHG4HitContainer.h \
  PHG4InEvent.h \
  PHG4Particle.h \
  PHG4Particlev1.h \
  PHG4Particlev2.h \
  PHG4ParticleGenerator.h \
  PHG4ParticleGeneratorBase.h \
  PHG4ParticleGeneratorVectorMeson.h \
  PHG4ParticleGun.h \
  PHG4PhenixDetector.h \
  PHG4PileupGenerator.h \
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
  PHG4VtxPointv1.h

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
	echo %_Dict_rdict.pcm

##############################################
# please add new classes in alphabetical order

PHG4_Dict.cc: \
  HepMCNodeReader.h \
  PHG4ConsistencyCheck.h \
  PHG4HeadReco.h \
  PHG4InEventCompress.h \
  PHG4InEventReadBack.h \
  PHG4InputFilter.h \
  PHG4PileupGenerator.h \
  PHG4SimpleEventGenerator.h \
  PHG4ParticleGun.h \
  PHG4ParticleGenerator.h \
  PHG4ParticleGeneratorBase.h \
  PHG4ParticleGeneratorVectorMeson.h \
  PHG4ParticleGeneratorD0.h \
  PHG4Reco.h \
  PHG4ScoringManager.h \
  PHG4Subsystem.h \
  PHG4TruthSubsystem.h \
  PHG4Utils.h \
  ReadEICFiles.h \
  PHG4LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(RINCLUDES) $^

PHG4AllHit_Dict.cc: \
  PHG4EventHeader.h \
  PHG4EventHeaderv1.h \
  PHG4HitReadBack.h \
  PHG4HitDefs.h \
  PHG4Hitv1.h \
  PHG4HitEval.h \
  PHG4HitContainer.h \
  PHG4InEvent.h \
  PHG4Shower.h \
  PHG4Showerv1.h \
  PHG4Particlev1.h \
  PHG4Particlev2.h \
  PHG4TruthInfoContainer.h \
  PHG4VtxPoint.h \
  PHG4VtxPointv1.h \
  PHG4AllHitLinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(RINCLUDES) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm 
