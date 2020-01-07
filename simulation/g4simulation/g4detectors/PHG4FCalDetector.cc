#include "PHG4FCalDetector.h"

#include "PHG4FCalSteppingAction.h"

#include <g4main/PHG4Detector.h>         // for PHG4Detector

#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Region.hh>  // for G4Region
#include <Geant4/G4String.hh>            // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Types.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cstdlib>   // for NULL, exit
#include <iostream>  // for stringstream, operator<<
#include <map>
#include <sstream>
#include <utility>  // for pair

class PHCompositeNode;

using namespace std;

PHG4FCalDetector::PHG4FCalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const string& name)
  : PHG4Detector(subsys, Node, name)
  , length(1.0 * m)
  , absorber_thickness(5.0 * mm)
  , scintillator_thickness(5.0 * cm)
  , nlayers(4)
  , segments_per_column(20)
  , segments_per_thickness(1)
  , z_position(100.0 * cm)
  , layer_separation(1.0 * mm)
  , AbsorberMaterial(nullptr)
  , ScintillatorMaterial(nullptr)
  , _region(nullptr)
{
}

G4Material* PHG4FCalDetector::SetMaterial(G4String material)
{
  // search the material by its name and assign if found
  if (G4Material* material_ptr = G4Material::GetMaterial(material))
  {
    return material_ptr;
  }
  else
    return 0;
}

unsigned int PHG4FCalDetector::computeIndex(unsigned int layer, G4double x, G4double y, G4double z, G4double& xcenter, G4double& ycenter, G4double& zcenter)
{
  double z_start = z_position + absorber_thickness * (((double) layer) + 1.) + layer_separation * (((double) layer) + 1.) + scintillator_thickness * (((double) layer));
  double z_width = scintillator_thickness / ((double) segments_per_thickness);
  unsigned int z_index = (unsigned int) ((z - z_start) / z_width);

  double xy_start = -0.5 * length;
  double xy_width = length / ((double) segments_per_column);

  unsigned int x_index = (unsigned int) ((x - xy_start) / xy_width);
  unsigned int y_index = (unsigned int) ((y - xy_start) / xy_width);

  zcenter = z_start + z_index * z_width + 0.5 * z_width;
  xcenter = xy_start + x_index * xy_width + 0.5 * xy_width;
  ycenter = xy_start + y_index * xy_width + 0.5 * xy_width;

  return (z_index * (segments_per_column * segments_per_column) + y_index * segments_per_column + x_index);
}

void PHG4FCalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  //   const G4MaterialTable* mattab = G4Material::GetMaterialTable();
  //   for(unsigned int i=0;i<mattab->size();i++)
  //   {
  //     cout<<mattab->at(i)->GetName()<<endl;
  //   }

  G4NistManager* man = G4NistManager::Instance();
  ScintillatorMaterial = man->FindOrBuildMaterial("G4_POLYSTYRENE");
  AbsorberMaterial = man->FindOrBuildMaterial("G4_Pb");

  if (!ScintillatorMaterial ||
      !AbsorberMaterial)
  {
    std::cout << "PHG4SvxDetector::Construct - Error: Can not set material" << std::endl;
    exit(-1);
  }

  G4VisAttributes* polystyreneVis = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  polystyreneVis->SetVisibility(true);
  polystyreneVis->SetForceSolid(true);

  G4VisAttributes* leadVis = new G4VisAttributes(G4Colour(0.925, 0.345, 0));
  leadVis->SetVisibility(true);
  leadVis->SetForceSolid(true);

  _region = new G4Region("FCALREGION");
  _region->SetRegionalSteppingAction(new PHG4FCalSteppingAction(this));

  G4double current_center_position = z_position;
  for (unsigned int layer = 0; layer < nlayers; layer++)
  {
    stringstream ss;

    ss.clear();
    ss.str("");
    ss << "FCalAbsSolid_" << layer;
    absorber_solid_[layer] = new G4Box(G4String(ss.str()), 0.5 * length, 0.5 * length, 0.5 * absorber_thickness);

    ss.clear();
    ss.str("");
    ss << "FCalAbsLogical_" << layer;
    absorber_logic_[layer] = new G4LogicalVolume(absorber_solid_[layer], AbsorberMaterial, G4String(ss.str()), 0, 0, 0);
    absorber_logic_[layer]->SetVisAttributes(leadVis);
    _region->AddRootLogicalVolume(absorber_logic_[layer]);

    ss.clear();
    ss.str("");
    ss << "FCalAbsPhysical_" << layer;

    absorber_physi_[layer] = new G4PVPlacement(0, G4ThreeVector(0. * cm, 0. * cm, current_center_position), absorber_logic_[layer], G4String(ss.str()), logicWorld, false, 0, false);

    current_center_position += 0.5 * absorber_thickness;
    current_center_position += 0.5 * scintillator_thickness;
    current_center_position += layer_separation;

    ss.clear();
    ss.str("");
    ss << "FCalScintSolid_" << layer;
    scintillator_solid_[layer] = new G4Box(G4String(ss.str()), 0.5 * length, 0.5 * length, 0.5 * scintillator_thickness);

    ss.clear();
    ss.str("");
    ss << "FCalScintLogical_" << layer;
    scintillator_logic_[layer] = new G4LogicalVolume(scintillator_solid_[layer], ScintillatorMaterial, G4String(ss.str()), 0, 0, 0);
    scintillator_logic_[layer]->SetVisAttributes(polystyreneVis);
    _region->AddRootLogicalVolume(scintillator_logic_[layer]);

    ss.clear();
    ss.str("");
    ss << "FCalScintPhysical_" << layer;

    scintillator_physi_[layer] = new G4PVPlacement(0, G4ThreeVector(0. * cm, 0. * cm, current_center_position), scintillator_logic_[layer], G4String(ss.str()), logicWorld, false, 0, false);

    current_center_position += 0.5 * absorber_thickness;
    current_center_position += 0.5 * scintillator_thickness;
    current_center_position += layer_separation;
  }
}

bool PHG4FCalDetector::isInScintillator(G4VPhysicalVolume* volume)
{
  //loop over the physical volumes and see if this is a match
  std::map<unsigned int, G4VPhysicalVolume*>::iterator vol_iter = scintillator_physi_.begin();
  for (; vol_iter != scintillator_physi_.end(); ++vol_iter)
  {
    if (vol_iter->second == volume)
      return true;
  }
  return false;
}

int PHG4FCalDetector::getScintillatorLayer(G4VPhysicalVolume* volume)
{
  //loop over the physical volumes and see if this is a match
  std::map<unsigned int, G4VPhysicalVolume*>::iterator vol_iter = scintillator_physi_.begin();
  for (; vol_iter != scintillator_physi_.end(); ++vol_iter)
  {
    if (vol_iter->second == volume)
      return vol_iter->first;
  }
  return -1;
}
