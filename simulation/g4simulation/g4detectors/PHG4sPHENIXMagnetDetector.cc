#include "PHG4sPHENIXMagnetDetector.h"

#include "PHG4sPHENIXMagnetDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4gdml/PHG4GDMLUtility.hh>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>

#include <TSystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4Polycone.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>

class G4Material;
class G4VSolid;
class PHCompositeNode;

//_______________________________________________________________________
PHG4sPHENIXMagnetDetector::PHG4sPHENIXMagnetDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam, const int detid)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4sPHENIXMagnetDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_GdmlConfig(PHG4GDMLUtility::GetOrMakeConfigNode(Node))
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_Layer(detid)
  , m_SuperDetector("NONE")
{
  assert(m_GdmlConfig);
}

//_______________________________________________________________________
int PHG4sPHENIXMagnetDetector::IsInsPHENIXMagnet(G4VPhysicalVolume* volume) const
{
  G4LogicalVolume* mylogvol = volume->GetLogicalVolume();

  if (m_LogicalVolSet.find(mylogvol) != m_LogicalVolSet.end())
  {
    return 1;
  }

  return 0;
}

//_______________________________________________________________________
void PHG4sPHENIXMagnetDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4sPHENIXMagnetDetector: Begin Construction" << std::endl;
  }

  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  G4Material* aluminum = GetDetectorMaterial("G4_Al");
  G4Material* copper = GetDetectorMaterial("G4_Cu");

  // Mother volume for solenoid

  G4Tubs* cryostat = SolenoidTubes(0);
  G4LogicalVolume* logCryostat =
      new G4LogicalVolume(cryostat, aluminum, "CRYOSTAT");
  GetDisplayAction()->AddVolume(logCryostat, "CRYOSTAT");
  m_LogicalVolSet.insert(logCryostat);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logCryostat, "CRYOSTAT", logicWorld, false, false, OverlapCheck());

  // Air (or vacuum?) inside cryostat

  G4Polycone* cryostatInterior = SolenoidPolycones(0);
  G4LogicalVolume* logCryostatInterior =
      new G4LogicalVolume(cryostatInterior, WorldMaterial, "CRYOINT");

  GetDisplayAction()->AddVolume(logCryostatInterior, "CRYOINT");
  m_LogicalVolSet.insert(logCryostatInterior);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logCryostatInterior, "CRYOINT", logCryostat, false, false, OverlapCheck());

  // Thermal shield

  G4Tubs* therm = SolenoidTubes(1);
  G4LogicalVolume* logTherm =
      new G4LogicalVolume(therm, aluminum, "THERM");
  GetDisplayAction()->AddVolume(logTherm, "THERM");
  m_LogicalVolSet.insert(logTherm);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logTherm, "THERM", logCryostatInterior, false, false, OverlapCheck());

  G4Tubs* thermvac = SolenoidTubes(2);
  G4LogicalVolume* logThermVac =
      new G4LogicalVolume(thermvac, WorldMaterial, "THERMVAC");
  GetDisplayAction()->AddVolume(logThermVac, "THERMVAC");
  m_LogicalVolSet.insert(logThermVac);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logThermVac, "THERMVAC", logTherm, false, false, OverlapCheck());

  // Coil support

  G4Polycone* coilSupport = SolenoidPolycones(1);
  G4LogicalVolume* logCoilSupport =
      new G4LogicalVolume(coilSupport, aluminum, "COILSUP");
  GetDisplayAction()->AddVolume(logCoilSupport, "COILSUP");
  m_LogicalVolSet.insert(logCoilSupport);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logCoilSupport, "COILSUP", logThermVac, false, false, OverlapCheck());

  // Coil

  G4Tubs* coil = SolenoidTubes(3);
  G4LogicalVolume* logCoil =
      new G4LogicalVolume(coil, aluminum, "COIL");
  GetDisplayAction()->AddVolume(logCoil, "COIL");
  m_LogicalVolSet.insert(logCoil);
  G4double zpos = 1.25 * cm;  // ifr_inert.dat file line 185
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, zpos),
                    logCoil, "COIL", logThermVac, false, false, OverlapCheck());

  // Tube connecting solenoid to cryotower

  G4Tubs* connector = CryoTubes(3);
  G4LogicalVolume* logConnector =
      new G4LogicalVolume(connector, aluminum, "CONNECTOR");
  GetDisplayAction()->AddVolume(logConnector, "CONNECTOR");
  m_LogicalVolSet.insert(logConnector);

  G4double ycon = 158.5 * cm;   // ifr_inert.dat file line 189
  G4double zcon = -197.5 * cm;  // ifr_inert.dat file line 189
  new G4PVPlacement(0, G4ThreeVector(0.0, ycon, zcon),
                    logConnector, "CONNECTOR", logicWorld, false, false, OverlapCheck());

  // Bus bar in connecting tube

  G4Box* busbarD = Block(14);
  G4LogicalVolume* logBusbarD =
      new G4LogicalVolume(busbarD, copper, "BusbarD", 0, 0, 0);
  GetDisplayAction()->AddVolume(logBusbarD, "BusbarD");
  m_LogicalVolSet.insert(logBusbarD);

  G4double yD = 158.5 * cm;   // ifr_inert.dat file line 75
  G4double zD = -197.5 * cm;  // ifr_inert.dat file line 75
  G4RotationMatrix rotationD;
  rotationD.rotateZ(90. * deg);
  rotationD.rotateX(90. * deg);
  G4Transform3D transformD(rotationD, G4ThreeVector(0.0, yD, zD));
  new G4PVPlacement(transformD,
                    logBusbarD, "BusbarD", logicWorld, false, false, OverlapCheck());

  return;
}

G4Tubs* PHG4sPHENIXMagnetDetector::SolenoidTubes(G4int itube) const
{
  G4double tubes_r_inner[8] = {142.0, 147.6, 148.1, 151.0, 12.44, 10.15, 0.0, 10.86};    // ifr_inert.dat file lines 181-189
  G4double tubes_r_outer[8] = {177.0, 168.0, 167.5, 155.08, 13.35, 10.50, 0.79, 12.37};  // ifr_inert.dat file lines 181-189
  G4double tubes_length[8] = {385.0, 371.0, 370.0, 351.29, 121.0, 121.0, 121.0, 9.98};   // ifr_inert.dat file lines 181-189

  G4double r_inner = tubes_r_inner[itube] * cm;
  G4double r_outer = tubes_r_outer[itube] * cm;
  G4double half_length = 0.5 * tubes_length[itube] * cm;

  return new G4Tubs("SolenoidTube", r_inner, r_outer, half_length,
                    0.0 * deg, 360.0 * deg);
}

G4Polycone* PHG4sPHENIXMagnetDetector::SolenoidPolycones(G4int ipolycone) const
{
  G4double zPlane[10];
  G4double rInner[10];
  G4double rOuter[10];

  G4int nzplanes = 10;  // ifr_inert.dat file lines 203-211
  // ifr_inert.dat file lines 203-211
  G4double polycones_zplane[2][10] = {{-188.5, -179.5, -172.0, -159.5, -150.2, 155.1, 164.4, 172.0, 179.5, 188.5},
                                      {-176.9, -174.4, -174.4, -164.0, -158.4, 163.4, 169.0, 176.9, 176.9, 181.9}};
  G4double polycones_r_inner[2][10] = {{145.5, 145.5, 143.0, 143.0, 143.0, 143.0, 143.0, 143.0, 145.5, 145.5},
                                       {151.0, 151.0, 155.1, 155.1, 155.1, 155.1, 155.1, 155.1, 151.0, 151.0}};
  G4double polycones_r_outer[2][10] = {{172.0, 172.0, 172.0, 172.0, 174.5, 174.5, 172.0, 172.0, 172.0, 172.0},
                                       {160.4, 160.4, 160.4, 160.4, 158.6, 158.6, 160.4, 160.4, 160.4, 160.4}};
  G4int iz;
  for (iz = 0; iz < nzplanes; iz++)
  {
    zPlane[iz] = polycones_zplane[ipolycone][iz] * cm;
    rInner[iz] = polycones_r_inner[ipolycone][iz] * cm;
    rOuter[iz] = polycones_r_outer[ipolycone][iz] * cm;
  }

  return new G4Polycone("SolenoidPolycone", 0.0 * deg, 360.0 * deg, nzplanes,
                        zPlane, rInner, rOuter);
}

G4Tubs* PHG4sPHENIXMagnetDetector::CryoTubes(G4int itube) const
{
  G4double tubes_r_inner[8] = {142.0, 147.6, 148.1, 151.0, 12.44, 10.15, 0.0, 10.86};    // ifr_inert.dat file lines 181-189
  G4double tubes_r_outer[8] = {177.0, 168.0, 167.5, 155.08, 13.35, 10.50, 0.79, 12.37};  // ifr_inert.dat file lines 181-189
  G4double tubes_length[8] = {385.0, 371.0, 370.0, 351.29, 121.0, 121.0, 121.0, 9.98};   // ifr_inert.dat file lines 181-189

  G4double r_inner = tubes_r_inner[itube + 4] * cm;
  G4double r_outer = tubes_r_outer[itube + 4] * cm;
  G4double half_length = 0.5 * tubes_length[itube + 4] * cm;

  return new G4Tubs("CryoTube", r_inner, r_outer, half_length,
                    0.0 * deg, 360.0 * deg);
}

G4Box* PHG4sPHENIXMagnetDetector::Block(G4int iblock) const
{
  G4double blocks_length[15] = {375.0, 138.9, 300.0, 315.0, 128.9, 128.9, 166.84, 160.84, 128.9, 128.9, 251.5, 121.0, 54.4, 14.6, 9.98};  // ifr_inert.dat file lines 60-75
  G4double blocks_width[15] = {172.52, 80.0, 290.0, 290.0, 104.9, 104.9, 104.9, 104.9, 114.9, 114.9, 15.0, 6.1, 6.1, 6.1, 6.1};           // ifr_inert.dat file lines 60-75
  G4double blocks_thickness[15] = {5.5, 15.0, 25.0, 25.0, 8.5, 8.5, 10.0, 10.0, 25.0, 25.0, 15.0, 3.4, 3.4, 3.4, 3.4};                    // ifr_inert.dat file lines 60-75

  G4double dx = 0.5 * blocks_length[iblock] * cm;
  G4double dy = 0.5 * blocks_width[iblock] * cm;
  G4double dz = 0.5 * blocks_thickness[iblock] * cm;

  return new G4Box("IfrBlock", dx, dy, dz);
}