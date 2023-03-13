AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libg4detectors_io.la \
  libg4detectors.la

AM_CXXFLAGS = `geant4-config --cflags`

AM_CPPFLAGS = \
  -DCGAL_DISABLE_ROUNDING_MATH_CHECK \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(XERCESROOT)/include \
  -I$(OPT_SPHENIX)/include

# set in configure.in to check if gcc version >= 4.8
#if GCC_GE_48
#  AM_CXXFLAGS = -std=c++11
#endif

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

libg4detectors_io_la_LIBADD = \
  -lphool \
  -lphparameter


#libgmp comes in via CGAL which is now header only
libg4detectors_la_LIBADD = \
  libg4detectors_io.la \
  -lffamodules \
  -lg4testbench \
  -lgmp \
  -lphg4gdml \
  -lphool  \
  -lSubsysReco

pkginclude_HEADERS = \
  BeamLineMagnetSubsystem.h \
  PHG4BbcSubsystem.h \
  PHG4BlockCellGeom.h \
  PHG4BlockCellGeomContainer.h \
  PHG4BlockGeom.h \
  PHG4BlockGeomContainer.h \
  PHG4BlockSubsystem.h \
  PHG4Cell.h \
  PHG4Cellv1.h \
  PHG4CellContainer.h \
  PHG4CellDefs.h \
  PHG4EventActionClearZeroEdep.h \
  PHG4ConeSubsystem.h \
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
  PHG4FullProjSpacalCellReco.h \
  PHG4GDMLSubsystem.h \
  PHG4HcalDefs.h \
  PHG4HcalCellReco.h \
  PHG4HcalSubsystem.h \
  PHG4InnerHcalSubsystem.h \
  PHG4OuterHcalSubsystem.h \
  PHG4PSTOFSubsystem.h \
  PHG4ScintillatorSlat.h \
  PHG4ScintillatorSlatv1.h \
  PHG4ScintillatorSlatContainer.h \
  PHG4ScintillatorSlatDefs.h \
  PHG4SectorConstructor.h \
  PHG4SectorSubsystem.h \
  PHG4SpacalSubsystem.h \
  PHG4StepStatusDecode.h \
  PHG4TpcCylinderGeom.h \
  PHG4TpcCylinderGeomContainer.h \
  PHG4ZDCDefs.h \
  PHG4ZDCSubsystem.h

ROOT_IO_DICTS = \
  PHG4BlockCellGeom_Dict.cc \
  PHG4BlockCellGeomContainer_Dict.cc \
  PHG4BlockGeom_Dict.cc \
  PHG4BlockGeomContainer_Dict.cc \
  PHG4BlockGeomv1_Dict.cc \
  PHG4Cell_Dict.cc \
  PHG4Cellv1_Dict.cc \
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
  PHG4ScintillatorSlatContainer_Dict.cc \
  PHG4TpcCylinderGeom_Dict.cc \
  PHG4TpcCylinderGeomContainer_Dict.cc 

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  PHG4BlockCellGeom_Dict_rdict.pcm \
  PHG4BlockCellGeomContainer_Dict_rdict.pcm \
  PHG4BlockGeom_Dict_rdict.pcm \
  PHG4BlockGeomContainer_Dict_rdict.pcm \
  PHG4BlockGeomv1_Dict_rdict.pcm \
  PHG4Cell_Dict_rdict.pcm \
  PHG4Cellv1_Dict_rdict.pcm \
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
  PHG4ScintillatorSlatContainer_Dict_rdict.pcm \
  PHG4TpcCylinderGeom_Dict_rdict.pcm \
  PHG4TpcCylinderGeomContainer_Dict_rdict.pcm 

libg4detectors_io_la_SOURCES = \
  $(ROOT_IO_DICTS) \
  PHG4BlockCellGeom.cc \
  PHG4BlockCellGeomContainer.cc \
  PHG4BlockGeom.cc \
  PHG4BlockGeomv1.cc \
  PHG4BlockGeomContainer.cc \
  PHG4Cell.cc \
  PHG4Cellv1.cc \
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
  PHG4ScintillatorSlatContainer.cc \
  PHG4TpcCylinderGeom.cc \
  PHG4TpcCylinderGeomContainer.cc 


libg4detectors_la_SOURCES = \
  BeamLineMagnetDetector.cc \
  BeamLineMagnetDisplayAction.cc \
  BeamLineMagnetSubsystem.cc \
  BeamLineMagnetSteppingAction.cc \
  PHG4GDMLDetector.cc \
  PHG4GDMLSubsystem.cc \
  PHG4BbcDetector.cc \
  PHG4BbcDisplayAction.cc \
  PHG4BbcSteppingAction.cc \
  PHG4BbcSubsystem.cc \
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
  PHG4ConeDisplayAction.cc \
  PHG4ConeSteppingAction.cc \
  PHG4ConeSubsystem.cc \
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
  PHG4FullProjSpacalDetector.cc \
  PHG4FullProjTiltedSpacalDetector.cc \
  PHG4FullProjSpacalCellReco.cc \
  PHG4GenHit.cc \
  PHG4HcalCellReco.cc \
  PHG4HcalDetector.cc \
  PHG4HcalSteppingAction.cc \
  PHG4HcalSubsystem.cc \
  PHG4InnerHcalDetector.cc \
  PHG4InnerHcalDisplayAction.cc \
  PHG4InnerHcalSteppingAction.cc \
  PHG4InnerHcalSubsystem.cc \
  PHG4OuterHcalDetector.cc \
  PHG4OuterHcalDisplayAction.cc \
  PHG4OuterHcalField.cc \
  PHG4OuterHcalFieldSetup.cc \
  PHG4OuterHcalSteppingAction.cc \
  PHG4OuterHcalSubsystem.cc \
  PHG4PSTOFDetector.cc \
  PHG4PSTOFSteppingAction.cc \
  PHG4PSTOFSubsystem.cc \
  PHG4SectorConstructor.cc \
  PHG4SectorDetector.cc \
  PHG4SectorDisplayAction.cc \
  PHG4SectorSteppingAction.cc \
  PHG4SectorSubsystem.cc \
  PHG4SpacalDetector.cc \
  PHG4SpacalDisplayAction.cc \
  PHG4SpacalSteppingAction.cc \
  PHG4SpacalSubsystem.cc \
  PHG4StepStatusDecode.cc \
  PHG4TrackStatusDecode.cc \
  PHG4ZDCDetector.cc \
  PHG4ZDCDisplayAction.cc \
  PHG4ZDCSteppingAction.cc \
  PHG4ZDCSubsystem.cc

bin_SCRIPTS = \
  CreateG4Subsystem.pl

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

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
