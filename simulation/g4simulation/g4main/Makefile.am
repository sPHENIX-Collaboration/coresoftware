AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags | sed 's/-I/-isystem/g'`

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
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(XERCESCROOT)/include \
  -isystem$(OPT_SPHENIX)/include

libphg4hit_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64


libphg4hit_la_LIBADD = \
  -lphool

libg4testbench_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `geant4-config --libs`

libg4testbench_la_LIBADD = \
  libphg4hit.la \
  -lboost_filesystem \
  -lffamodules \
  -lfun4all \
  -lg4decayer \
  -lgsl \
  -lgslcblas \
  -lHepMC \
  -lphfield \
  -lphgeom \
  -lphg4gdml \
  -lphhepmc \
  -lphparameter

# I/O dictionaries have to exist for root5 and root6. For ROOT6 we need
# pcm files in addition. If someone can figure out how to make a list
# so this list of dictionaries is transformed into a list of pcm files
# following a simple naming convention, please change this accordingly and
# let me know
ROOTHITDICTS = \
  EicEventHeader_Dict.cc \
  EicEventHeaderv1_Dict.cc \
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
  PHG4VtxPointv2_Dict.cc \
  PHG4TruthInfoContainer_Dict.cc

# this is a tweak to install files in $(libdir), automake refuses
# to install other files in libdir, this construct gets around this
pcmdir = $(libdir)
nobase_dist_pcm_DATA = $(ROOTHITDICTS:.cc=_rdict.pcm)

libphg4hit_la_SOURCES = \
  $(ROOTHITDICTS) \
  EicEventHeader.cc \
  EicEventHeaderv1.cc \
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
  PHG4Shower.cc \
  PHG4Showerv1.cc \
  PHG4TruthInfoContainer.cc \
  PHG4VtxPoint.cc \
  PHG4VtxPointv1.cc \
  PHG4VtxPointv2.cc

libg4testbench_la_SOURCES = \
  Fun4AllMessenger.cc \
  G4TBMagneticFieldSetup.cc \
  G4TBFieldMessenger.cc \
  Fun4AllDstPileupInputManager.cc \
  Fun4AllDstPileupMerger.cc \
  Fun4AllSingleDstPileupInputManager.cc \
  HepMCNodeReader.cc \
  PHG4ConsistencyCheck.cc \
  PHG4DisplayAction.cc \
  PHG4Detector.cc \
  PHG4HitReadBack.cc \
  PHG4EtaParameterization.cc \
  PHG4EtaPhiParameterization.cc \
  PHG4HeadReco.cc \
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
  PHG4PhenixStackingAction.cc \
  PHG4PhenixSteppingAction.cc \
  PHG4PhenixTrackingAction.cc \
  PHG4PileupGenerator.cc \
  PHG4PrimaryGeneratorAction.cc \
  PHG4ProcessMap.cc \
  PHG4ProcessMapPhysics.cc \
  PHG4Reco.cc \
  PHG4RegionInformation.cc \
  PHG4ScoringManager.cc \
  PHG4SimpleEventGenerator.cc \
  PHG4StackingAction.cc \
  PHG4SteppingAction.cc \
  PHG4Subsystem.cc \
  PHG4TrackUserInfoV1.cc \
  PHG4TruthEventAction.cc \
  PHG4TruthSubsystem.cc \
  PHG4TruthTrackingAction.cc \
  PHG4UIsession.cc \
  PHG4Utils.cc \
  PHG4VertexSelection.cc \
  CosmicSpray.cc


##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  EicEventHeader.h \
  Fun4AllDstPileupInputManager.h \
  Fun4AllDstPileupMerger.h \
  Fun4AllSingleDstPileupInputManager.h \
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
  PHG4MCProcessDefs.h \
  PHG4Particle.h \
  PHG4Particlev1.h \
  PHG4Particlev2.h \
  PHG4Particlev3.h \
  PHG4MagneticField.h \
  PHG4ParticleGenerator.h \
  PHG4ParticleGeneratorBase.h \
  PHG4ParticleGeneratorD0.h \
  PHG4ParticleGeneratorVectorMeson.h \
  PHG4ParticleGun.h \
  PHG4PhenixDetector.h \
  PHG4PileupGenerator.h \
  PHG4PrimaryGeneratorAction.h \
  PHG4ProcessMap.h \
  PHG4ProcessMapPhysics.h \
  PHG4Reco.h \
  PHG4RegionInformation.h \
  PHG4ScoringManager.h \
  PHG4SimpleEventGenerator.h \
  PHG4Shower.h \
  PHG4Showerv1.h \
  PHG4StackingAction.h \
  PHG4SteppingAction.h \
  PHG4Subsystem.h \
  PHG4TrackingAction.h \
  PHG4TrackUserInfoV1.h \
  PHG4TruthInfoContainer.h \
  PHG4TruthSubsystem.h \
  PHG4Units.h \
  PHG4Utils.h \
  PHG4VertexSelection.h \
  PHG4VtxPoint.h \
  PHG4VtxPointv1.h \
  PHG4VtxPointv2.h \
  ReadEICFiles.h \
  CosmicSpray.h \
  EcoMug.h

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
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
