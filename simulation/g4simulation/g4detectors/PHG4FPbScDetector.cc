#include "PHG4FPbScDetector.h"

#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Region.hh>         // for G4Region
#include <Geant4/G4SystemOfUnits.hh>  // for cm
#include <Geant4/G4ThreeVector.hh>    // for G4ThreeVector
#include <Geant4/G4Types.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cstdlib>   // for NULL, exit
#include <iostream>  // for stringstream, operator<<
#include <map>
#include <sstream>
#include <utility>  // for pair

class PHCompositeNode;

using namespace std;

PHG4FPbScDetector::PHG4FPbScDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& nam)
  : PHG4Detector(subsys, Node, nam)
  , tower_cross_section(5.535 * cm)
  , segments_per_column(12 * 6)
  , segments_per_height(12 * 3)
  , length(tower_cross_section * segments_per_column)
  , height(tower_cross_section * segments_per_height)
  , absorber_thickness(0.15 * cm)
  , scintillator_thickness(0.4 * cm)
  , nlayers(66)
  , x_position(0.0 * cm)
  , y_position(0.0 * cm)
  , z_position(400.0 * cm)
  , layer_separation(0.0)
  , AbsorberMaterial(nullptr)
  , ScintillatorMaterial(nullptr)
  , _region(nullptr)
{
}

G4Material* PHG4FPbScDetector::SetMaterial(G4String material)
{
  // search the material by its name and assign if found
  if (G4Material* material_ptr = G4Material::GetMaterial(material))
  {
    return material_ptr;
  }
  else
    return 0;
}

unsigned int PHG4FPbScDetector::computeIndex(unsigned int layer, G4double x, G4double y, G4double z, G4double& xcenter, G4double& ycenter, G4double& zcenter)
{
  double z_start = z_position + absorber_thickness * (((double) layer) + 1.) + layer_separation * (((double) layer) + 1.) + scintillator_thickness * (((double) layer));
  double z_width = scintillator_thickness;
  unsigned int z_index = (unsigned int) ((z - z_start) / z_width);

  double x_start = x_position - 0.5 * length;
  double y_start = y_position - 0.5 * height;
  double xy_width = tower_cross_section;

  unsigned int x_index = (unsigned int) ((x - x_start) / xy_width);
  unsigned int y_index = (unsigned int) ((y - y_start) / xy_width);

  zcenter = z;  //z_start + z_index*z_width + 0.5*z_width;
  xcenter = x;  //x_start + x_index*xy_width + 0.5*xy_width;
  ycenter = y;  //y_start + y_index*xy_width + 0.5*xy_width;

  return (z_index * (segments_per_column * segments_per_height) + y_index * segments_per_column + x_index);
}

void PHG4FPbScDetector::ConstructMe(G4LogicalVolume* logicWorld)
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
    std::cout << "PHG4FPbScDetector::Construct - Error: Can not set material" << std::endl;
    exit(-1);
  }

  G4VisAttributes* polystyreneVis = new G4VisAttributes(G4Colour::Blue());
  polystyreneVis->SetVisibility(true);
  polystyreneVis->SetForceSolid(true);

  G4VisAttributes* leadVis = new G4VisAttributes(G4Colour(0.925, 0.345, 0));
  leadVis->SetVisibility(true);
  leadVis->SetForceSolid(true);

  _region = new G4Region("FPBSCREGION");
  //  _region->SetRegionalSteppingAction(new PHG4FPbScSteppingAction(this));

  G4double current_center_position = z_position;
  for (unsigned int layer = 0; layer < nlayers; layer++)
  {
    stringstream ss;

    ss.clear();
    ss.str("");
    ss << "FPbScAbsSolid_" << layer;
    absorber_solid_[layer] = new G4Box(G4String(ss.str()),
                                       0.5 * length,
                                       0.5 * height,
                                       0.5 * absorber_thickness);

    ss.clear();
    ss.str("");
    ss << "FPbScAbsLogical_" << layer;
    absorber_logic_[layer] = new G4LogicalVolume(absorber_solid_[layer],
                                                 AbsorberMaterial,
                                                 G4String(ss.str()), 0, 0, 0);
    absorber_logic_[layer]->SetVisAttributes(leadVis);
    //    _region->AddRootLogicalVolume(absorber_logic_[layer]);

    ss.clear();
    ss.str("");
    ss << "FPbScAbsPhysical_" << layer;

    absorber_physi_[layer] = new G4PVPlacement(0,
                                               G4ThreeVector(x_position, y_position,
                                                             current_center_position),
                                               absorber_logic_[layer],
                                               G4String(ss.str()),
                                               logicWorld, false, 0, false);

    current_center_position += 0.5 * absorber_thickness;
    current_center_position += 0.5 * scintillator_thickness;

    ss.clear();
    ss.str("");
    ss << "FPbScScintSolid_" << layer;
    scintillator_solid_[layer] = new G4Box(G4String(ss.str()),
                                           0.5 * length, 0.5 * height,
                                           0.5 * scintillator_thickness);

    ss.clear();
    ss.str("");
    ss << "FPbScScintLogical_" << layer;
    scintillator_logic_[layer] = new G4LogicalVolume(scintillator_solid_[layer],
                                                     ScintillatorMaterial, G4String(ss.str()), 0, 0, 0);
    scintillator_logic_[layer]->SetVisAttributes(polystyreneVis);
    //    _region->AddRootLogicalVolume(scintillator_logic_[layer]);

    ss.clear();
    ss.str("");
    ss << "FPbScScintPhysical_" << layer;

    scintillator_physi_[layer] = new G4PVPlacement(0, G4ThreeVector(x_position, y_position, current_center_position),
                                                   scintillator_logic_[layer],
                                                   G4String(ss.str()),
                                                   logicWorld, false, 0, false);

    current_center_position += 0.5 * scintillator_thickness;
    current_center_position += layer_separation;
    current_center_position += 0.5 * absorber_thickness;
  }
}

bool PHG4FPbScDetector::isInScintillator(G4VPhysicalVolume* volume)
{
  //loop over the physical volumes and see if this is a match
  std::map<unsigned int, G4VPhysicalVolume*>::iterator vol_iter = scintillator_physi_.begin();
  for (; vol_iter != scintillator_physi_.end(); ++vol_iter)
  {
    if (vol_iter->second == volume)
    {
      return true;
    }
  }
  return false;
}

int PHG4FPbScDetector::getScintillatorLayer(G4VPhysicalVolume* volume)
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
