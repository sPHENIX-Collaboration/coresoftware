AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libg4detectors_io.la \
    libg4detectors.la

AM_CPPFLAGS = \
    -DCGAL_DISABLE_ROUNDING_MATH_CHECK \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(ROOTSYS)/include \
    -I$(G4_MAIN)/include \
    -I$(OPT_SPHENIX)/include

# set in configure.in to check if gcc version >= 4.8
#if GCC_GE_48
#  AM_CXXFLAGS = -std=c++11
#endif

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libg4detectors_io_la_LIBADD = \
  -lphool \
  -lphparameter


libg4detectors_la_LIBADD = \
  libg4detectors_io.la \
  -lphg4gdml \
  -lphool  \
  -lCGAL \
  -lSubsysReco \
  -lg4testbench 

pkginclude_HEADERS = \
  ePHENIXRICHConstruction.h \
  PHG4BlockCellGeom.h \
  PHG4BlockCellGeomContainer.h \
  PHG4BlockGeom.h \
  PHG4BlockGeomContainer.h \
  PHG4BlockSubsystem.h \
  PHG4Cell.h \
  PHG4Cellv1.h \
  PHG4Cellv2.h \
  PHG4CellContainer.h \
  PHG4CellDefs.h \
  PHG4EventActionClearZeroEdep.h \
  PHG4ConeSubsystem.h \
  PHG4CrystalCalorimeterSubsystem.h \
  PHG4CylinderCell.h \
  PHG4CylinderCellv1.h \
  PHG4CylinderCellv2.h \
  PHG4CylinderCellv3.h \
  PHG4CylinderCellContainer.h \
  PHG4CylinderGeom.h \
  PHG4CylinderGeomv1.h \
  PHG4CylinderGeomv2.h \
  PHG4CylinderGeomv3.h \
  PHG4CylinderGeomv4.h \
  PHG4CylinderGeom_Spacalv1.h \
  PHG4CylinderGeom_Spacalv2.h \
  PHG4CylinderGeom_Spacalv3.h \
  PHG4CylinderGeomContainer.h \
  PHG4CylinderCellDefs.h \
  PHG4CylinderCellGeom.h \
  PHG4CylinderCellGeom_Spacalv1.h \
  PHG4CylinderCellGeomContainer.h \
  PHG4CylinderCellReco.h \
  PHG4CylinderSubsystem.h \
  PHG4DetectorSubsystem.h \
  PHG4DetectorGroupSubsystem.h \
  PHG4ForwardCalCellReco.h \
  PHG4ForwardEcalSubsystem.h \
  PHG4ForwardHcalSubsystem.h \
  PHG4FullProjSpacalCellReco.h \
  PHG4HcalDefs.h \
  PHG4HcalCellReco.h \
  PHG4InnerHcalSubsystem.h \
  PHG4OuterHcalSubsystem.h \
  PHG4PrototypeHcalDefs.h \
  PHG4Prototype2HcalCellReco.h \
  PHG4Prototype2OuterHcalSubsystem.h \
  PHG4Prototype2InnerHcalSubsystem.h \
  PHG4Prototype3InnerHcalSubsystem.h \
  PHG4PSTOFSubsystem.h \
  PHG4RICHSubsystem.h \
  PHG4ScintillatorSlat.h \
  PHG4ScintillatorSlatContainer.h \
  PHG4ScintillatorSlatDefs.h \
  PHG4SectorConstructor.h \
  PHG4SectorSubsystem.h \
  PHG4SpacalSubsystem.h \
  PHG4SpacalPrototypeSubsystem.h \
  PHG4StepStatusDecode.h 

# I/O dictionaries have to exist for root5 and root6. For ROOT6 we need
# pcm files in addition. If someone can figure out how to make a list
# so this list of dictionaries is transformed into a list of pcm files
# following a simple naming convention, please change this accordingly and
# let me know
ROOT_IO_DICTS = \
  PHG4BlockCellGeom_Dict.cc \
  PHG4BlockCellGeomContainer_Dict.cc \
  PHG4BlockGeom_Dict.cc \
  PHG4BlockGeomContainer_Dict.cc \
  PHG4BlockGeomv1_Dict.cc \
  PHG4Cell_Dict.cc \
  PHG4Cellv1_Dict.cc \
  PHG4Cellv2_Dict.cc \
  PHG4CellContainer_Dict.cc \
  PHG4CylinderGeom_Dict.cc \
  PHG4CylinderGeomv1_Dict.cc \
  PHG4CylinderGeomv2_Dict.cc \
  PHG4CylinderGeomv3_Dict.cc \
  PHG4CylinderGeomv4_Dict.cc \
  PHG4CylinderGeom_Spacalv1_Dict.cc \
  PHG4CylinderGeom_Spacalv2_Dict.cc \
  PHG4CylinderGeom_Spacalv3_Dict.cc \
  PHG4CylinderGeomContainer_Dict.cc \
  PHG4CylinderCell_Dict.cc \
  PHG4CylinderCellv1_Dict.cc \
  PHG4CylinderCellv2_Dict.cc \
  PHG4CylinderCellv3_Dict.cc \
  PHG4CylinderCellContainer_Dict.cc \
  PHG4CylinderCellGeom_Dict.cc \
  PHG4CylinderCellGeom_Spacalv1_Dict.cc \
  PHG4CylinderCellGeomContainer_Dict.cc \
  PHG4ScintillatorSlat_Dict.cc \
  PHG4ScintillatorSlatv1_Dict.cc \
  PHG4ScintillatorSlatContainer_Dict.cc
# for root6 we need pcm and dictionaries but only for
# i/o classes. For root5 we need only dictionaries but
# those for i/o and classes available on the cmd line
# MAKEROOT6 is set in the configure.ac, true for root6
if MAKEROOT6
# this is a tweak to install files in $(libdir), automake refuses
# to install other files in libdir, this construct gets around this
pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  PHG4BlockCellGeom_Dict_rdict.pcm \
  PHG4BlockCellGeomContainer_Dict_rdict.pcm \
  PHG4BlockGeom_Dict_rdict.pcm \
  PHG4BlockGeomContainer_Dict_rdict.pcm \
  PHG4BlockGeomv1_Dict_rdict.pcm \
  PHG4Cell_Dict_rdict.pcm \
  PHG4Cellv1_Dict_rdict.pcm \
  PHG4Cellv2_Dict_rdict.pcm \
  PHG4CellContainer_Dict_rdict.pcm \
  PHG4CylinderGeom_Dict_rdict.pcm \
  PHG4CylinderGeomv1_Dict_rdict.pcm \
  PHG4CylinderGeomv2_Dict_rdict.pcm \
  PHG4CylinderGeomv3_Dict_rdict.pcm \
  PHG4CylinderGeomv4_Dict_rdict.pcm \
  PHG4CylinderGeom_Spacalv1_Dict_rdict.pcm \
  PHG4CylinderGeom_Spacalv2_Dict_rdict.pcm \
  PHG4CylinderGeom_Spacalv3_Dict_rdict.pcm \
  PHG4CylinderGeomContainer_Dict_rdict.pcm \
  PHG4CylinderCell_Dict_rdict.pcm \
  PHG4CylinderCellv1_Dict_rdict.pcm \
  PHG4CylinderCellv2_Dict_rdict.pcm \
  PHG4CylinderCellv3_Dict_rdict.pcm \
  PHG4CylinderCellContainer_Dict_rdict.pcm \
  PHG4CylinderCellGeom_Dict_rdict.pcm \
  PHG4CylinderCellGeom_Spacalv1_Dict_rdict.pcm \
  PHG4CylinderCellGeomContainer_Dict_rdict.pcm \
  PHG4ScintillatorSlat_Dict_rdict.pcm \
  PHG4ScintillatorSlatv1_Dict_rdict.pcm \
  PHG4ScintillatorSlatContainer_Dict_rdict.pcm
else
  ROOT5_IO_DICTS = \
  PHG4CellDefs_Dict.cc \
  PHG4ScintillatorSlatDefs_Dict.cc
  ROOT5_DICTS = \
  ePHENIXRICHConstruction_Dict.cc \
  PHG4GDMLSubsystem_Dict.cc \
  PHG4BeamlineMagnetSubsystem_Dict.cc \
  PHG4BlockCellReco_Dict.cc \
  PHG4BlockSubsystem_Dict.cc \
  PHG4CEmcTestBeamSubsystem_Dict.cc \
  PHG4ConeSubsystem_Dict.cc \
  PHG4CrystalCalorimeterSubsystem_Dict.cc \
  PHG4CylinderSubsystem_Dict.cc \
  PHG4CylinderCellReco_Dict.cc \
  PHG4DetectorGroupSubsystem_Dict.cc \
  PHG4DetectorSubsystem_Dict.cc \
  PHG4EnvelopeSubsystem_Dict.cc \
  PHG4FCalSubsystem_Dict.cc \
  PHG4ForwardEcalSubsystem_Dict.cc \
  PHG4ForwardHcalSubsystem_Dict.cc \
  PHG4ForwardCalCellReco_Dict.cc \
  PHG4FPbScSubsystem_Dict.cc \
  PHG4GenHit_Dict.cc \
  PHG4HcalCellReco_Dict.cc \
  PHG4HcalSubsystem_Dict.cc \
  PHG4HcalPrototypeSubsystem_Dict.cc \
  PHG4InnerHcalSubsystem_Dict.cc \
  PHG4OuterHcalSubsystem_Dict.cc \
  PHG4Prototype2InnerHcalSubsystem_Dict.cc \
  PHG4Prototype2HcalCellReco_Dict.cc \
  PHG4Prototype2OuterHcalSubsystem_Dict.cc \
  PHG4Prototype3InnerHcalSubsystem_Dict.cc \
  PHG4RICHSubsystem_Dict.cc \
  PHG4SectorConstructor_Dict.cc \
  PHG4SectorSubsystem_Dict.cc \
  PHG4PSTOFSubsystem_Dict.cc \
  PHG4FullProjSpacalCellReco_Dict.cc \
  PHG4SpacalSubsystem_Dict.cc \
  PHG4SpacalPrototypeSubsystem_Dict.cc \
  PHG4mRICHSubsystem_Dict.cc

endif

libg4detectors_io_la_SOURCES = \
  $(ROOT_IO_DICTS) \
  $(ROOT5_IO_DICTS) \
  PHG4BlockCellGeom.cc \
  PHG4BlockCellGeomContainer.cc \
  PHG4BlockGeom.cc \
  PHG4BlockGeomv1.cc \
  PHG4BlockGeomContainer.cc \
  PHG4Cell.cc \
  PHG4Cellv1.cc \
  PHG4Cellv2.cc \
  PHG4CellContainer.cc \
  PHG4CellDefs.cc \
  PHG4CylinderGeom.cc \
  PHG4CylinderGeomv1.cc \
  PHG4CylinderGeomv2.cc \
  PHG4CylinderGeomv3.cc \
  PHG4CylinderGeomv4.cc \
  PHG4CylinderGeom_Spacalv1.cc \
  PHG4CylinderGeom_Spacalv2.cc \
  PHG4CylinderGeom_Spacalv3.cc \
  PHG4CylinderGeom_Spacalv3_2015.cc \
  PHG4CylinderGeomContainer.cc \
  PHG4CylinderCellv1.cc \
  PHG4CylinderCellv2.cc \
  PHG4CylinderCellv3.cc \
  PHG4CylinderCellContainer.cc \
  PHG4CylinderCellGeom.cc \
  PHG4CylinderCellGeom_Spacalv1.cc \
  PHG4CylinderCellGeomContainer.cc \
  PHG4ScintillatorSlatDefs.cc \
  PHG4ScintillatorSlatv1.cc \
  PHG4ScintillatorSlatContainer.cc

libg4detectors_la_SOURCES = \
  $(ROOT5_DICTS) \
  PHG4GDMLDetector.cc \
  PHG4GDMLSubsystem.cc \
  PHG4BeamlineMagnetDetector.cc \
  PHG4BeamlineMagnetSubsystem.cc \
  PHG4BlockCellReco.cc \
  PHG4BlockDetector.cc \
  PHG4BlockDisplayAction.cc \
  PHG4BlockSteppingAction.cc \
  PHG4BlockSubsystem.cc \
  PHG4CEmcTestBeamDetector.cc \
  PHG4CEmcTestBeamSteppingAction.cc \
  PHG4CEmcTestBeamSubsystem.cc \
  PHG4ConeDetector.cc \
  PHG4ConeRegionSteppingAction.cc \
  PHG4ConeSteppingAction.cc \
  PHG4ConeSubsystem.cc \
  PHG4CrystalCalorimeterDetector.cc \
  PHG4CrystalCalorimeterSteppingAction.cc \
  PHG4CrystalCalorimeterSubsystem.cc \
  PHG4CylinderDetector.cc \
  PHG4CylinderDisplayAction.cc \
  PHG4CylinderSubsystem.cc \
  PHG4CylinderCellReco.cc \
  PHG4CylinderSteppingAction.cc \
  PHG4DetectorGroupSubsystem.cc \
  PHG4DetectorSubsystem.cc \
  PHG4EnvelopeDetector.cc \
  PHG4EnvelopeSteppingAction.cc \
  PHG4EnvelopeSubsystem.cc \
  PHG4EventActionClearZeroEdep.cc \
  PHG4FCalDetector.cc \
  PHG4FCalSteppingAction.cc \
  PHG4FCalSubsystem.cc \
  PHG4ForwardEcalDetector.cc \
  PHG4EICForwardEcalDetector.cc \
  PHG4ForwardEcalSteppingAction.cc \
  PHG4ForwardEcalSubsystem.cc \
  PHG4ForwardHcalDetector.cc \
  PHG4ForwardHcalSteppingAction.cc \
  PHG4ForwardHcalSubsystem.cc \
  PHG4ForwardCalCellReco.cc \
  PHG4FPbScDetector.cc \
  PHG4FPbScSteppingAction.cc \
  PHG4FPbScRegionSteppingAction.cc \
  PHG4FPbScSubsystem.cc \
  PHG4GenHit.cc \
  PHG4HcalCellReco.cc \
  PHG4HcalDetector.cc \
  PHG4HcalSteppingAction.cc \
  PHG4HcalSubsystem.cc \
  PHG4HcalPrototypeDetector.cc \
  PHG4HcalPrototypeDetectorMessenger.cc \
  PHG4HcalPrototypeSteppingAction.cc \
  PHG4HcalPrototypeSubsystem.cc \
  PHG4InnerHcalDetector.cc \
  PHG4InnerHcalSteppingAction.cc \
  PHG4InnerHcalSubsystem.cc \
  PHG4OuterHcalDetector.cc \
  PHG4OuterHcalField.cc \
  PHG4OuterHcalFieldSetup.cc \
  PHG4OuterHcalSteppingAction.cc \
  PHG4OuterHcalSubsystem.cc \
  PHG4ProjCrystalCalorimeterDetector.cc \
  PHG4Prototype2InnerHcalDetector.cc \
  PHG4Prototype2InnerHcalSteppingAction.cc \
  PHG4Prototype2InnerHcalSubsystem.cc \
  PHG4Prototype2HcalCellReco.cc \
  PHG4Prototype2OuterHcalDetector.cc \
  PHG4Prototype2OuterHcalSteppingAction.cc \
  PHG4Prototype2OuterHcalSubsystem.cc \
  PHG4Prototype3InnerHcalDetector.cc \
  PHG4Prototype3InnerHcalSteppingAction.cc \
  PHG4Prototype3InnerHcalSubsystem.cc \
  PHG4RICHDetector.cc \
  PHG4RICHSteppingAction.cc \
  PHG4RICHSubsystem.cc \
  ePHENIXRICHConstruction.cc \
  PHG4SectorConstructor.cc \
  PHG4SectorDetector.cc \
  PHG4SectorSteppingAction.cc \
  PHG4SectorSubsystem.cc \
  PHG4StepStatusDecode.cc \
  PHG4TrackStatusDecode.cc \
  PHG4PSTOFDetector.cc \
  PHG4PSTOFSteppingAction.cc \
  PHG4PSTOFSubsystem.cc \
  PHG4SpacalDetector.cc \
  PHG4FullProjSpacalDetector.cc \
  PHG4FullProjTiltedSpacalDetector.cc \
  PHG4FullProjSpacalCellReco.cc \
  PHG4SpacalPrototypeDetector.cc \
  PHG4SpacalSteppingAction.cc \
  PHG4SpacalPrototypeSteppingAction.cc \
  PHG4SpacalSubsystem.cc \
  PHG4SpacalPrototypeSubsystem.cc \
  PHG4mRICHDetector.cc \
  PHG4mRICHSteppingAction.cc \
  PHG4mRICHSubsystem.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_g4detectors \
  testexternals_g4detectors_io

testexternals_g4detectors_SOURCES = testexternals.cc
testexternals_g4detectors_LDADD = libg4detectors.la

testexternals_g4detectors_io_SOURCES = testexternals.cc
testexternals_g4detectors_io_LDADD = libg4detectors_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
