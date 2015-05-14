AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libg4detectors.la

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

INCLUDES = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(ROOTSYS)/include \
    -I${G4_MAIN}/include \
    -I${G4_MAIN}/include/Geant4 \
    -I/direct/phenix+u/nfeege/sphenixsw/devel/install/g4main/include

libg4detectors_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lg4testbench \
  -lCGAL

pkginclude_HEADERS = \
  PHG4BlockGeom.h \
  PHG4BlockGeomContainer.h \
  PHG4CylinderCell.h \
  PHG4CylinderCellv1.h \
  PHG4CylinderCellv2.h \
  PHG4CylinderCellContainer.h \
  PHG4CylinderGeom.h \
  PHG4CylinderGeomContainer.h \
  PHG4CylinderCellDefs.h \
  PHG4CylinderCellGeom.h \
  PHG4CylinderCellGeomContainer.h

libg4detectors_la_SOURCES = \
  PHG4AccDetector.cc \
  PHG4AccSteppingAction.cc \
  PHG4AccSubsystem.cc \
  PHG4BlockCellGeom.cc \
  PHG4BlockCellGeomContainer.cc \
  PHG4BlockCellReco.cc \
  PHG4BlockDetector.cc \
  PHG4BlockEventAction.cc \
  PHG4BlockGeom.cc \
  PHG4BlockGeomv1.cc \
  PHG4BlockGeomContainer.cc \
  PHG4BlockRegionSteppingAction.cc \
  PHG4BlockSteppingAction.cc \
  PHG4BlockSubsystem.cc \
  PHG4CEmcTestBeamDetector.cc \
  PHG4CEmcTestBeamSteppingAction.cc \
  PHG4CEmcTestBeamSubsystem.cc \
  PHG4EventActionClearZeroEdep.cc \
  PHG4ConeDetector.cc \
  PHG4ConeEventAction.cc \
  PHG4ConeRegionSteppingAction.cc \
  PHG4ConeSteppingAction.cc \
  PHG4ConeSubsystem.cc \
  PHG4CrystalCalorimeterDetector.cc \
  PHG4CrystalCalorimeterSteppingAction.cc \
  PHG4CrystalCalorimeterSubsystem.cc \
  PHG4CylinderDetector.cc \
  PHG4CylinderEventAction.cc \
  PHG4CylinderGeom.cc \
  PHG4CylinderGeomv1.cc \
  PHG4CylinderGeomv2.cc \
  PHG4CylinderGeomv3.cc \
  PHG4CylinderGeomv4.cc \
  PHG4CylinderGeom_Spacalv1.cc \
  PHG4CylinderGeom_Spacalv2.cc \
  PHG4CylinderGeomContainer.cc \
  PHG4CylinderRegionSteppingAction.cc \
  PHG4CylinderSubsystem.cc \
  PHG4CylinderCell.cc \
  PHG4CylinderCellv1.cc \
  PHG4CylinderCellv2.cc \
  PHG4CylinderCellContainer.cc \
  PHG4CylinderCellGeom.cc \
  PHG4CylinderCellGeomContainer.cc \
  PHG4CylinderCellReco.cc \
  PHG4CylinderSteppingAction.cc \
  PHG4Detector_Dict.cc \
  PHG4FCalDetector.cc \
  PHG4FCalSteppingAction.cc \
  PHG4FCalSubsystem.cc \
  PHG4FPbScDetector.cc \
  PHG4FPbScEventAction.cc \
  PHG4FPbScSteppingAction.cc \
  PHG4FPbScRegionSteppingAction.cc \
  PHG4FPbScSubsystem.cc \
  PHG4GenHit.cc \
  PHG4HcalCellReco.cc \
  PHG4HcalDetector.cc \
  PHG4HcalSteppingAction.cc \
  PHG4HcalSubsystem.cc \
  PHG4HcalTestBeamDetector.cc \
  PHG4HcalTestBeamSteppingAction.cc \
  PHG4HcalTestBeamSubsystem.cc \
  PHG4HcalPrototypeDetector.cc \
  PHG4HcalPrototypeDetectorMessenger.cc \
  PHG4HcalPrototypeSteppingAction.cc \
  PHG4HcalPrototypeSubsystem.cc \
  PHG4InnerHcalDetector.cc \
  PHG4InnerHcalSteppingAction.cc \
  PHG4InnerHcalSubsystem.cc \
  PHG4SlatCellReco.cc \
  PHG4OuterHcalDetector.cc \
  PHG4OuterHcalSteppingAction.cc \
  PHG4OuterHcalSubsystem.cc \
  PHG4RICHDetector.cc \
  PHG4RICHSteppingAction.cc \
  PHG4RICHSubsystem.cc \
  ePHENIXRICHConstruction.cc \
  PHG4SectorConstructor.cc \
  PHG4SectorDetector.cc \
  PHG4SectorEventAction.cc \
  PHG4SectorSteppingAction.cc \
  PHG4SectorSubsystem.cc \
  PHG4SiliconTrackerCellReco.cc \
  PHG4SiliconTrackerDetector.cc \
  PHG4SiliconTrackerSteppingAction.cc \
  PHG4SiliconTrackerSubsystem.cc \
  PHG4SpacalDetector.cc \
  PHG4ProjSpacalDetector.cc \
  PHG4SpacalSteppingAction.cc \
  PHG4SpacalSubsystem.cc

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

PHG4Detector_Dict.cc : \
  ePHENIXRICHConstruction.h \
  PHG4AccSubsystem.h \
  PHG4BlockCellGeom.h \
  PHG4BlockCellGeomContainer.h \
  PHG4BlockCellReco.h \
  PHG4BlockSubsystem.h \
  PHG4BlockGeom.h \
  PHG4BlockGeomv1.h \
  PHG4BlockGeomContainer.h \
  PHG4CEmcTestBeamSubsystem.h \
  PHG4ConeSubsystem.h \
  PHG4CrystalCalorimeterSubsystem.h \
  PHG4CylinderCell.h \
  PHG4CylinderCellv1.h \
  PHG4CylinderCellv2.h \
  PHG4CylinderCellContainer.h \
  PHG4CylinderCellGeom.h \
  PHG4CylinderCellGeomContainer.h \
  PHG4CylinderCellReco.h \
  PHG4CylinderGeom.h \
  PHG4CylinderGeomv1.h \
  PHG4CylinderGeomv2.h \
  PHG4CylinderGeomv3.h \
  PHG4CylinderGeomv4.h \
  PHG4CylinderGeomContainer.h \
  PHG4CylinderGeom_Spacalv1.h \
  PHG4CylinderGeom_Spacalv2.h \
  PHG4CylinderSubsystem.h \
  PHG4FCalSubsystem.h \
  PHG4FPbScSubsystem.h \
  PHG4GenHit.h \
  PHG4HcalCellReco.h \
  PHG4HcalSubsystem.h \
  PHG4HcalTestBeamSubsystem.h \
  PHG4HcalPrototypeSubsystem.h \
  PHG4InnerHcalSubsystem.h \
  PHG4OuterHcalSubsystem.h \
  PHG4RICHSubsystem.h \
  PHG4SectorConstructor.h \
  PHG4SectorSubsystem.h \
  PHG4SiliconTrackerSubsystem.h \
  PHG4SiliconTrackerCellReco.h \
  PHG4SlatCellReco.h \
  PHG4SpacalSubsystem.h \
  PHG4DetectorLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

