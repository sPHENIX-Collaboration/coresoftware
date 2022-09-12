// $Id: $

/*!
 * \file PHG4GDMLDetector.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4GDMLDetector.h"

#include <g4main/PHG4Detector.h>  // for PHG4Detector
#include <g4main/PHG4Subsystem.h>
#include <g4main/PHG4Utils.h>

#include <phparameter/PHParameters.h>

#include <g4gdml/PHG4GDMLConfig.hh>
#include <g4gdml/PHG4GDMLUtility.hh>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4VisAttributes.hh>

// Xerces has shadowed variables
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4GDMLReadStructure.hh>  // for G4GDMLReadStructure
#pragma GCC diagnostic pop

#include <CLHEP/Units/SystemOfUnits.h>  // for cm, degree

#include <cassert>
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, basic_ostream
#include <memory>
#include <vector>  // for vector, vector<>::iterator

PHG4GDMLDetector::PHG4GDMLDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam, PHParameters* parameters)
  : PHG4Detector(subsys, Node, dnam)
  , m_GDMPath(parameters->get_string_param("GDMPath"))
  , m_TopVolName(parameters->get_string_param("TopVolName"))
  , m_placeX(parameters->get_double_param("place_x") * cm)
  , m_placeY(parameters->get_double_param("place_y") * cm)
  , m_placeZ(parameters->get_double_param("place_z") * cm)
  , m_rotationX(parameters->get_double_param("rot_x") * degree)
  , m_rotationY(parameters->get_double_param("rot_y") * degree)
  , m_rotationZ(parameters->get_double_param("rot_z") * degree)
  , m_skipDSTGeometryExport(parameters->get_int_param("skip_DST_geometry_export"))
{
  if (m_skipDSTGeometryExport)
  {
    gdml_config = PHG4GDMLUtility::GetOrMakeConfigNode(Node);
    assert(gdml_config);
  }
}

PHG4GDMLDetector::~PHG4GDMLDetector() = default;

void

PHG4GDMLDetector::Print(const std::string& /*what*/) const
{
  std::cout << "PHG4GDMLDetector::" << GetName() << " - import " << m_TopVolName << " from " << m_GDMPath << " with shift "
            << m_placeX << ","
            << m_placeY << ","
            << m_placeZ << "cm and rotation "
            << m_rotationX << ","
            << m_rotationY << ","
            << m_rotationZ << "rad" << std::endl;
}

void PHG4GDMLDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    std::cout << " PHG4MapsDetector::Construct:";
    Print();
    //      std::cout << std::endl;
  }

  //===================================
  // Import the stave physical volume here
  //===================================

  // import the staves from the gemetry file
  std::unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.SetOverlapCheck(OverlapCheck());
  //  gdmlParser.Read(m_GDMPath, false);
  gdmlParser.Read(m_GDMPath, OverlapCheck());

  //  G4AssemblyVolume* av_ITSUStave = reader->GetAssembly(assemblyname);

  G4LogicalVolume* vol = reader->GetVolume(m_TopVolName);

  if (not vol)
  {
    std::cout << "PHG4GDMLDetector::Construct - Fatal Error - failed to find G4LogicalVolume " << m_TopVolName << " - Print: ";
    Print();
    exit(121);
  }
  PHG4Subsystem* mysys = GetMySubsystem();
  mysys->SetLogicalVolume(vol);

  G4RotationMatrix* rotm = new G4RotationMatrix();
  rotm->rotateX(m_rotationX);
  rotm->rotateY(m_rotationY);
  rotm->rotateZ(m_rotationZ);
  G4ThreeVector placeVec(m_placeX, m_placeY, m_placeZ);

  //  av_ITSUStave->MakeImprint(trackerenvelope, Tr, 0, OverlapCheck());

  G4PVPlacement* gdml_phys =
      new G4PVPlacement(rotm, placeVec,
                        vol,
                        G4String(GetName()),
                        logicWorld, false, 0, OverlapCheck());
  SetDisplayProperty(vol);

  if (m_skipDSTGeometryExport)
  {
    assert(gdml_config);
    gdml_config->exclude_physical_vol(gdml_phys);
  }
}

void PHG4GDMLDetector::SetDisplayProperty(G4AssemblyVolume* av)
{
  //  std::cout <<"SetDisplayProperty - G4AssemblyVolume w/ TotalImprintedVolumes "<<av->TotalImprintedVolumes()
  //   <<"/"<<av->GetImprintsCount()<<std::endl;

  std::vector<G4VPhysicalVolume*>::iterator it = av->GetVolumesIterator();

  int nDaughters = av->TotalImprintedVolumes();
  for (int i = 0; i < nDaughters; ++i, ++it)
  {
    //  std::cout <<"SetDisplayProperty - AV["<<i<<"] = "<<(*it)->GetName()<<std::endl;
    G4VPhysicalVolume* pv = (*it);

    G4LogicalVolume* worldLogical = pv->GetLogicalVolume();
    SetDisplayProperty(worldLogical);
  }
}

void PHG4GDMLDetector::SetDisplayProperty(G4LogicalVolume* lv)
{
  std::string material_name(
      lv->GetMaterial()->GetName());

  if (Verbosity() >= 5)
    std::cout << "SetDisplayProperty - LV " << lv->GetName() << " built with "
              << material_name << std::endl;

  G4VisAttributes* matVis = new G4VisAttributes();
  if (material_name.find("SI") != std::string::npos)
  {
    PHG4Utils::SetColour(matVis, "G4_Si");
    matVis->SetVisibility(true);
    matVis->SetForceSolid(true);
    if (Verbosity() >= 5)
      std::cout << "SetDisplayProperty - LV " << lv->GetName() << " display with G4_Si" << std::endl;
  }
  else if (material_name.find("KAPTON") != std::string::npos)
  {
    PHG4Utils::SetColour(matVis, "G4_KAPTON");
    matVis->SetVisibility(true);
    matVis->SetForceSolid(true);
    if (Verbosity() >= 5)
      std::cout << "SetDisplayProperty - LV " << lv->GetName() << " display with G4_KAPTON" << std::endl;
  }
  else if (material_name.find("ALUMINUM") != std::string::npos)
  {
    PHG4Utils::SetColour(matVis, "G4_Al");
    matVis->SetVisibility(true);
    matVis->SetForceSolid(true);
    if (Verbosity() >= 5)
      std::cout << "SetDisplayProperty - LV " << lv->GetName() << " display with G4_Al" << std::endl;
  }
  else if (material_name.find("Carbon") != std::string::npos)
  {
    matVis->SetColour(0.5, 0.5, 0.5, .25);
    matVis->SetVisibility(true);
    matVis->SetForceSolid(true);
    if (Verbosity() >= 5)
      std::cout << "SetDisplayProperty - LV " << lv->GetName() << " display with Gray" << std::endl;
  }
  else if (material_name.find("M60J3K") != std::string::npos)
  {
    matVis->SetColour(0.25, 0.25, 0.25, .25);
    matVis->SetVisibility(true);
    matVis->SetForceSolid(true);
    if (Verbosity() >= 5)
      std::cout << "SetDisplayProperty - LV " << lv->GetName() << " display with Gray" << std::endl;
  }
  else if (material_name.find("WATER") != std::string::npos)
  {
    matVis->SetColour(0.0, 0.5, 0.0, .25);
    matVis->SetVisibility(true);
    matVis->SetForceSolid(true);
    if (Verbosity() >= 5)
      std::cout << "SetDisplayProperty - LV " << lv->GetName() << " display with WATER" << std::endl;
  }
  else
  {
    matVis->SetColour(.2, .2, .7, .25);
    matVis->SetVisibility(true);
    matVis->SetForceSolid(true);
  }
  lv->SetVisAttributes(matVis);

  int nDaughters = lv->GetNoDaughters();
  for (int i = 0; i < nDaughters; ++i)
  {
    G4VPhysicalVolume* pv = lv->GetDaughter(i);

    // std::cout <<"SetDisplayProperty - PV["<<i<<"] = "<<pv->GetName()<<std::endl;

    G4LogicalVolume* worldLogical = pv->GetLogicalVolume();
    SetDisplayProperty(worldLogical);
  }
}
