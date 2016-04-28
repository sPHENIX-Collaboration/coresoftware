AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libg4detectors_io.la \
    libg4detectors.la

AM_CPPFLAGS = -DCGAL_DISABLE_ROUNDING_MATH_CHECK

# set in configure.in to check if gcc version >= 4.8
if GCC_GE_48
  AM_CXXFLAGS = -std=c++11
endif

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

INCLUDES = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(ROOTSYS)/include \
    -I${G4_MAIN}/include \
    -I${G4_MAIN}/include/Geant4

libg4detectors_io_la_LIBADD = \
  -lphool \
  -lboost_filesystem \
  -lboost_system \
  -lpdbcalBase

libg4detectors_la_LIBADD = \
  -lphool  \
  libg4detectors_io.la \
  -lCGAL \
  -lSubsysReco \
  -lg4testbench 

pkginclude_HEADERS = \
  PHG4Parameters.h \
  PHG4BlockGeom.h \
  PHG4BlockGeomContainer.h \
  PHG4CylinderCell.h \
  PHG4CylinderCellv1.h \
  PHG4CylinderCell_Spacalv1.h \
  PHG4CylinderCellv2.h \
  PHG4CylinderCellv3.h \
  PHG4CylinderCellContainer.h \
  PHG4CylinderGeom.h \
  PHG4CylinderGeom_Spacalv1.h \
  PHG4CylinderGeom_Spacalv2.h \
  PHG4CylinderGeom_Spacalv3.h \
  PHG4CylinderGeomContainer.h \
  PHG4CylinderCellDefs.h \
  PHG4CylinderCellGeom.h \
  PHG4CylinderCellGeom_Spacalv1.h \
  PHG4CylinderCellGeomContainer.h \
  PHG4ScintillatorSlat.h \
  PHG4ScintillatorSlatContainer.h \
  PHG4ScintillatorSlatDefs.h

libg4detectors_io_la_SOURCES = \
  PHG4Parameters.cc \
  PHG4Parameters_Dict.cc \
  PHG4BlockCellGeom.cc \
  PHG4BlockCellGeom_Dict.cc \
  PHG4BlockCellGeomContainer.cc \
  PHG4BlockCellGeomContainer_Dict.cc \
  PHG4BlockGeom.cc \
  PHG4BlockGeom_Dict.cc \
  PHG4BlockGeomv1.cc \
  PHG4BlockGeomv1_Dict.cc \
  PHG4BlockGeomContainer.cc \
  PHG4BlockGeomContainer_Dict.cc \
  PHG4CylinderGeom.cc \
  PHG4CylinderGeom_Dict.cc \
  PHG4CylinderGeomv1.cc \
  PHG4CylinderGeomv1_Dict.cc \
  PHG4CylinderGeomv2.cc \
  PHG4CylinderGeomv2_Dict.cc \
  PHG4CylinderGeomv3.cc \
  PHG4CylinderGeomv3_Dict.cc \
  PHG4CylinderGeomv4.cc \
  PHG4CylinderGeomv4_Dict.cc \
  PHG4CylinderGeom_Spacalv1.cc \
  PHG4CylinderGeom_Spacalv1_Dict.cc \
  PHG4CylinderGeom_Spacalv2.cc \
  PHG4CylinderGeom_Spacalv2_Dict.cc \
  PHG4CylinderGeom_Spacalv3.cc \
  PHG4CylinderGeom_Spacalv3_Dict.cc \
  PHG4CylinderGeomContainer.cc \
  PHG4CylinderGeomContainer_Dict.cc \
  PHG4CylinderCell.cc \
  PHG4CylinderCell_Dict.cc \
  PHG4CylinderCellv1.cc \
  PHG4CylinderCellv1_Dict.cc \
  PHG4CylinderCell_Spacalv1.cc \
  PHG4CylinderCell_Spacalv1_Dict.cc \
  PHG4CylinderCellv2.cc \
  PHG4CylinderCellv2_Dict.cc \
  PHG4CylinderCellv3.cc \
  PHG4CylinderCellv3_Dict.cc \
  PHG4CylinderCellContainer.cc \
  PHG4CylinderCellContainer_Dict.cc \
  PHG4CylinderCellGeom.cc \
  PHG4CylinderCellGeom_Dict.cc \
  PHG4CylinderCellGeom_Spacalv1.cc \
  PHG4CylinderCellGeom_Spacalv1_Dict.cc \
  PHG4CylinderCellGeomContainer.cc \
  PHG4CylinderCellGeomContainer_Dict.cc \
  PHG4ScintillatorSlat.cc \
  PHG4ScintillatorSlat_Dict.cc \
  PHG4ScintillatorSlatContainer.cc \
  PHG4ScintillatorSlatContainer_Dict.cc

libg4detectors_la_SOURCES = \
  PHG4BlockCellReco.cc \
  PHG4BlockCellReco_Dict.cc \
  PHG4BlockDetector.cc \
  PHG4BlockRegionSteppingAction.cc \
  PHG4BlockSteppingAction.cc \
  PHG4BlockSubsystem.cc \
  PHG4BlockSubsystem_Dict.cc \
  PHG4CEmcTestBeamDetector.cc \
  PHG4CEmcTestBeamSteppingAction.cc \
  PHG4CEmcTestBeamSubsystem.cc \
  PHG4CEmcTestBeamSubsystem_Dict.cc \
  PHG4EventActionClearZeroEdep.cc \
  PHG4ConeDetector.cc \
  PHG4ConeRegionSteppingAction.cc \
  PHG4ConeSteppingAction.cc \
  PHG4ConeSubsystem.cc \
  PHG4ConeSubsystem_Dict.cc \
  PHG4CrystalCalorimeterDetector.cc \
  PHG4CrystalCalorimeterSteppingAction.cc \
  PHG4CrystalCalorimeterSubsystem.cc \
  PHG4CrystalCalorimeterSubsystem_Dict.cc \
  PHG4CylinderDetector.cc \
  PHG4CylinderRegionSteppingAction.cc \
  PHG4CylinderSubsystem.cc \
  PHG4CylinderSubsystem_Dict.cc \
  PHG4CylinderCellReco.cc \
  PHG4CylinderCellReco_Dict.cc \
  PHG4CylinderSteppingAction.cc \
  PHG4EnvelopeDetector.cc \
  PHG4EnvelopeSteppingAction.cc \
  PHG4EnvelopeSubsystem.cc \
  PHG4EnvelopeSubsystem_Dict.cc \
  PHG4FCalDetector.cc \
  PHG4FCalSteppingAction.cc \
  PHG4FCalSubsystem.cc \
  PHG4FCalSubsystem_Dict.cc \
  PHG4ForwardEcalDetector.cc \
  PHG4EICForwardEcalDetector.cc \
  PHG4ForwardEcalSteppingAction.cc \
  PHG4ForwardEcalSubsystem.cc \
  PHG4ForwardEcalSubsystem_Dict.cc \
  PHG4ForwardHcalDetector.cc \
  PHG4ForwardHcalSteppingAction.cc \
  PHG4ForwardHcalSubsystem.cc \
  PHG4ForwardHcalSubsystem_Dict.cc \
  PHG4ForwardCalCellReco.cc \
  PHG4ForwardCalCellReco_Dict.cc \
  PHG4FPbScDetector.cc \
  PHG4FPbScSteppingAction.cc \
  PHG4FPbScRegionSteppingAction.cc \
  PHG4FPbScSubsystem.cc \
  PHG4FPbScSubsystem_Dict.cc \
  PHG4GenHit.cc \
  PHG4GenHit_Dict.cc \
  PHG4HcalCellReco.cc \
  PHG4HcalCellReco_Dict.cc \
  PHG4HcalDetector.cc \
  PHG4HcalSteppingAction.cc \
  PHG4HcalSubsystem.cc \
  PHG4HcalSubsystem_Dict.cc \
  PHG4HcalPrototypeDetector.cc \
  PHG4HcalPrototypeDetectorMessenger.cc \
  PHG4HcalPrototypeSteppingAction.cc \
  PHG4HcalPrototypeSubsystem.cc \
  PHG4HcalPrototypeSubsystem_Dict.cc \
  PHG4HcalPrototype2Detector.cc \
  PHG4HcalPrototype2DetectorMessenger.cc \
  PHG4HcalPrototype2SteppingAction.cc \
  PHG4HcalPrototype2Subsystem.cc \
  PHG4HcalPrototype2Subsystem_Dict.cc \
  PHG4InnerHcalDetector.cc \
  PHG4InnerHcalSteppingAction.cc \
  PHG4InnerHcalSubsystem.cc \
  PHG4InnerHcal_Dict.C \
  PHG4ScintillatorSlat.cc \
  PHG4ScintillatorSlat_Dict.cc \
  PHG4ScintillatorSlatDefs.cc \
  PHG4ScintillatorSlatDefs_Dict.cc \
  PHG4ScintillatorSlatv1.cc \
  PHG4ScintillatorSlatv1_Dict.cc \
  PHG4SlatCellReco.cc \
  PHG4SlatCellReco_Dict.cc \
  PHG4OuterHcalDetector.cc \
  PHG4OuterHcalField.cc \
  PHG4OuterHcalFieldSetup.cc \
  PHG4OuterHcalSteppingAction.cc \
  PHG4OuterHcalSubsystem.cc \
  PHG4OuterHcalSubsystem_Dict.cc \
  PHG4ProjCrystalCalorimeterDetector.cc \
  PHG4Prototype2InnerHcalDetector.cc \
  PHG4Prototype2InnerHcalSteppingAction.cc \
  PHG4Prototype2InnerHcalSubsystem.cc \
  PHG4Prototype2InnerHcalSubsystem_Dict.C \
  PHG4Prototype2HcalCellReco.cc \
  PHG4Prototype2HcalCellReco_Dict.cc \
  PHG4Prototype2OuterHcalDetector.cc \
  PHG4Prototype2OuterHcalSteppingAction.cc \
  PHG4Prototype2OuterHcalSubsystem.cc \
  PHG4Prototype2OuterHcalSubsystem_Dict.C \
  PHG4RICHDetector.cc \
  PHG4RICHSteppingAction.cc \
  PHG4RICHSubsystem.cc \
  PHG4RICHSubsystem_Dict.cc \
  ePHENIXRICHConstruction.cc \
  ePHENIXRICHConstruction_Dict.cc \
  PHG4SectorConstructor.cc \
  PHG4SectorConstructor_Dict.cc \
  PHG4SectorDetector.cc \
  PHG4SectorSteppingAction.cc \
  PHG4SectorSubsystem.cc \
  PHG4SectorSubsystem_Dict.cc \
  PHG4SiliconTrackerCellReco.cc \
  PHG4SiliconTrackerCellReco_Dict.cc \
  PHG4SiliconTrackerDetector.cc \
  PHG4SiliconTrackerSteppingAction.cc \
  PHG4SiliconTrackerSubsystem.cc \
  PHG4SiliconTrackerSubsystem_Dict.cc \
  PHG4SpacalDetector.cc \
  PHG4ProjSpacalDetector.cc \
  PHG4FullProjSpacalDetector.cc \
  PHG4FullProjSpacalCellReco.cc \
  PHG4FullProjSpacalCellReco_Dict.cc \
  PHG4SpacalPrototypeDetector.cc \
  PHG4SpacalSteppingAction.cc \
  PHG4SpacalPrototypeSteppingAction.cc \
  PHG4SpacalSubsystem.cc \
  PHG4SpacalSubsystem_Dict.cc \
  PHG4SpacalPrototypeSubsystem.cc \
  PHG4SpacalPrototypeSubsystem_Dict.cc \
  PHG4CylinderCellTPCReco.cc \
  PHG4CylinderCellTPCReco_Dict.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

PHG4InnerHcal_Dict.C: \
  PHG4InnerHcalSubsystem.h \
  PHG4InnerHcalLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4detectors

testexternals_g4detectors_SOURCES = testexternals.cc
testexternals_g4detectors_LDADD = libg4detectors.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
