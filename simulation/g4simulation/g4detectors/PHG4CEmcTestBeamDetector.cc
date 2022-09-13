#include "PHG4CEmcTestBeamDetector.h"

#include <g4main/PHG4Detector.h>  // for PHG4Detector

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>  // for G4double
#include <Geant4/G4VisAttributes.hh>

#include <algorithm>  // for copy
#include <cmath>      // for cos, sin, NAN, acos, atan
#include <iostream>   // for operator<<, ostringstream
#include <sstream>

class G4Material;
class G4VSolid;
class PHCompositeNode;

using namespace std;

static double no_overlap = 0.0001 * cm;  // added safety margin against overlaps by using same boundary between volumes

PHG4CEmcTestBeamDetector::PHG4CEmcTestBeamDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , gap(0.25 * mm)
  , place_in_x(0 * cm)
  , place_in_y(0 * cm)
  , place_in_z(0 * cm)
  , plate_x(135 * mm)
  , plate_z(135 * mm)
  , x_rot(0)
  , y_rot(0)
  , z_rot(0)
  , alpha(NAN)
  , inner_radius(NAN)
  , outer_radius(NAN)
  , tower_angular_coverage(NAN)
  , cemc_angular_coverage(NAN)
  , active_scinti_fraction(0.78)
  , sandwiches_per_tower(12)
  ,  // 12 tungsten/scintillator fiber snadwiches per tower
  num_towers(7)
  , active(0)
  , absorberactive(0)
  , layer(lyr)
  , blackhole(0)
{
  w_dimension[0] = plate_x;
  w_dimension[1] = 0.5 * mm;
  w_dimension[2] = plate_z;
  sc_dimension[0] = plate_x;
  sc_dimension[1] = 1 * mm;
  sc_dimension[2] = plate_z;
  sandwich_thickness = 2 * w_dimension[1] + sc_dimension[1];  // two tungsten plats, one scintillator
  for (int i = 0; i < 4; i++)
  {
    sandwich_vol.push_back(nullptr);
  }
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4CEmcTestBeamDetector::IsInCEmcTestBeam(G4VPhysicalVolume* volume) const
{
  if (active && volume == sandwich_vol[2])
  {
    return 1;
  }
  if (absorberactive && (volume == sandwich_vol[0] || volume == sandwich_vol[1]))
  {
    return -1;
  }
  if (absorberactive && sandwich_vol[3] && volume == sandwich_vol[3])
  {
    return -1;
  }
  return 0;
}

void PHG4CEmcTestBeamDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  CalculateGeometry();
  recoConsts* rc = recoConsts::instance();
  G4Material* Air = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  G4VSolid* cemc_tub = new G4Tubs("CEmcTub", inner_radius - 2 * no_overlap, outer_radius + 2 * no_overlap, (w_dimension[2] + 2 * no_overlap) / 2., 0, cemc_angular_coverage);
  G4LogicalVolume* cemc_log = new G4LogicalVolume(cemc_tub, Air, G4String("CEmc"), nullptr, nullptr, nullptr);

  G4RotationMatrix cemc_rotm;
  // put our cemc at center displacement in x
  double radius_at_center = inner_radius + (outer_radius - inner_radius) / 2.;
  double xcenter = radius_at_center * cos(cemc_angular_coverage / 2.);
  double ycenter = radius_at_center * sin(cemc_angular_coverage / 2.);
  cemc_rotm.rotateX(x_rot);
  cemc_rotm.rotateY(y_rot);
  cemc_rotm.rotateZ(z_rot);
  new G4PVPlacement(G4Transform3D(cemc_rotm, G4ThreeVector(place_in_x - xcenter, place_in_y - ycenter, place_in_z)), cemc_log, "CEmc", logicWorld, false, false, OverlapCheck());
  G4VSolid* tower_tub = new G4Tubs("TowerTub", inner_radius - no_overlap, outer_radius + no_overlap, (w_dimension[2] + no_overlap) / 2., 0, tower_angular_coverage);
  G4LogicalVolume* tower_log = new G4LogicalVolume(tower_tub, Air, G4String("CEmcTower"), nullptr, nullptr, nullptr);
  //  G4VisAttributes* towerVisAtt = new G4VisAttributes();
  // towerVisAtt->SetVisibility(true);
  // towerVisAtt->SetForceSolid(true);
  // towerVisAtt->SetColour(G4Colour::Blue());
  //  tower_log->SetVisAttributes(towerVisAtt);

  ostringstream tower_vol_name;
  for (int i = 0; i < 7; i++)
  {
    tower_vol_name << "CEmcTower_" << i;
    double phi = -i * tower_angular_coverage;
    G4RotationMatrix* tower_rotm = new G4RotationMatrix();
    tower_rotm->rotateZ(phi * rad);
    new G4PVPlacement(tower_rotm, G4ThreeVector(0, 0, 0), tower_log, tower_vol_name.str(), cemc_log, false, i, OverlapCheck());
    tower_vol_name.str("");
  }
  ConstructTowerVolume(tower_log);
  return;
}

// here we build up the tower from the sandwichs (tungsten + scintillator)
int PHG4CEmcTestBeamDetector::ConstructTowerVolume(G4LogicalVolume* tower_log)
{
  recoConsts* rc = recoConsts::instance();
  G4Material* Air = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  G4VSolid* sandwich_box = new G4Box("Sandwich_box",
                                     w_dimension[0] / 2., sandwich_thickness / 2., w_dimension[2] / 2.);

  G4LogicalVolume* sandwich_log = new G4LogicalVolume(sandwich_box,
                                                      Air,
                                                      G4String("CEmcSandwich"),
                                                      nullptr, nullptr, nullptr);
  G4VisAttributes* sandwichVisAtt = new G4VisAttributes();
  sandwichVisAtt->SetVisibility(true);
  sandwichVisAtt->SetForceSolid(true);
  sandwichVisAtt->SetColour(G4Colour::White());
  sandwich_log->SetVisAttributes(sandwichVisAtt);
  ostringstream sandwich_name;
  for (int i = 0; i < 12; i++)
  {
    G4RotationMatrix* sandwich_rotm = new G4RotationMatrix();
    double phi = -i * alpha;
    sandwich_rotm->rotateZ(phi * rad);
    sandwich_name << "CEmcSandwich_" << i;
    double xshift = cos(phi) * (inner_radius + (outer_radius - inner_radius) / 2.);
    double yshift = -sin(phi) * (inner_radius + (outer_radius - inner_radius) / 2.);
    // we need to shift everything up by sandwich_thickness/2, calculate the shift in x
    // if we push the sandwich up by sandwich_thickness/2
    double xcorr = atan(phi) * sandwich_thickness / 2.;
    double ycorr = sandwich_thickness / 2.;

    new G4PVPlacement(sandwich_rotm, G4ThreeVector(xshift + xcorr, yshift + ycorr, 0),
                      sandwich_log,
                      sandwich_name.str(),
                      tower_log, false, i, OverlapCheck());
    sandwich_name.str("");
  }
  ConstructSandwichVolume(sandwich_log);  // put W and scinti into sandwich
  return 0;
}

// here we put a single tungsten + scintillator sandwich together
int PHG4CEmcTestBeamDetector::ConstructSandwichVolume(G4LogicalVolume* sandwich)
{
  vector<G4LogicalVolume*> block_logic;
  G4Material* AbsorberMaterial = GetDetectorMaterial("G4_W");
  G4Material* ScintiMaterial = GetDetectorMaterial("G4_POLYSTYRENE");

  if (active_scinti_fraction > 1 || active_scinti_fraction < 0)
  {
    cout << "invalid active scintillator fraction " << active_scinti_fraction
         << " try between 0 and 1" << endl;
  }

  double sc_active_thickness = sc_dimension[1] * active_scinti_fraction;
  double sc_passive_thickness = sc_dimension[1] - sc_active_thickness;

  G4VSolid* block_w = new G4Box("Tungsten_box",
                                w_dimension[0] / 2., w_dimension[1] / 2., w_dimension[2] / 2.);
  G4VSolid* block_sc = new G4Box("Scinti_box",
                                 sc_dimension[0] / 2., sc_active_thickness / 2., sc_dimension[2] / 2.);
  G4VSolid* block_passive_sc = nullptr;
  block_logic.push_back(new G4LogicalVolume(block_w,
                                            AbsorberMaterial,
                                            "Plate_log_W",
                                            nullptr, nullptr, nullptr));
  block_logic.push_back(new G4LogicalVolume(block_sc,
                                            ScintiMaterial,
                                            "Plate_log_Sc",
                                            nullptr, nullptr, nullptr));
  G4VisAttributes* matVis = new G4VisAttributes();
  G4VisAttributes* matVis1 = new G4VisAttributes();
  matVis->SetVisibility(true);
  matVis->SetForceSolid(true);
  matVis->SetColour(G4Colour::Red());
  matVis1->SetVisibility(true);
  matVis1->SetForceSolid(true);
  matVis1->SetColour(G4Colour::Green());
  block_logic[0]->SetVisAttributes(matVis);
  block_logic[1]->SetVisAttributes(matVis1);

  // here is the math to set the positions of those tiles inside the sandwich. The center of the coordinate system is in the center
  // of the sandwich (w -> thickness of tungsten layer, sc_a -> thickness of actvie scinti, sc_p -> thickness of passive scinti)
  // distance to the bottom of the sandwich:
  // -(w + sc_a + sc_p + w)/2 = -w/2 -sc_a/2 -sc_p/2 -w/2

  // implement the absorber at the bottom and at the top of the sandwich
  // moving up by w/2 for the lowest layer tungsten:
  // -w/2 -sc_a/2 -sc_p/2 -w/2 + w/2 = -w/2 -sc_a/2 -sc_p/2 = -(w+sc)/2
  // where sc_a+sc_p = sc = scintillator thickness (sc_dimension[1])

  sandwich_vol[0] = new G4PVPlacement(nullptr, G4ThreeVector(0, -(w_dimension[1] + sc_dimension[1]) / 2., 0),
                                      block_logic[0],
                                      "CEmc_W_plate_down",
                                      sandwich, false, 0, OverlapCheck());

  // top of the sandwich - starting at the bottom and add the tungsten and scintillator layer and half tungsten
  // -w/2 -sc/2 -w/2 + w + sc + w/2 = +sc/2 +w/2 = (w+sc)/2
  sandwich_vol[1] = new G4PVPlacement(nullptr, G4ThreeVector(0, (w_dimension[1] + sc_dimension[1]) / 2., 0),
                                      block_logic[0],
                                      "CEmc_W_plate_up",
                                      sandwich, false, 0, OverlapCheck());

  // since we split the scintillator into an active and passive part to accomodate
  // for the scintillator fibers not occupying 100% of the volume.
  // this is mocked up by 2 separate layers of plastic,  active scintillator is the middle layer
  // again going to the bottom of the sandwich box (sc_a -> active sc, sc_p -> passive sc)
  //  -w/2 -sc_a/2 -sc_p/2 -w/2
  // -(w + sc_a + sc_p + w)/2, adding the w layer and half of the sc_a to get to the middle of the sc_a layer:
  // -w/2 -sc_a/2 - sc_p/2 -w/2 + w + sc_a/2 = -sc_p/2
  // if fraction of active sc is 1, sc_passive_thickness is zero
  sandwich_vol[2] = new G4PVPlacement(nullptr, G4ThreeVector(0, -sc_passive_thickness / 2., 0),
                                      block_logic[1],
                                      "CEmc_Scinti_plate",
                                      sandwich, false, 0, OverlapCheck());

  if (sc_passive_thickness > 0)
  {
    G4VisAttributes* matVis2 = new G4VisAttributes();
    matVis1->SetVisibility(true);
    matVis1->SetForceSolid(true);
    matVis1->SetColour(G4Colour::Blue());
    block_passive_sc = new G4Box("Passive_Scinti_box",
                                 sc_dimension[0] / 2., sc_passive_thickness / 2., sc_dimension[2] / 2.);
    block_logic.push_back(new G4LogicalVolume(block_passive_sc,
                                              ScintiMaterial,
                                              "Plate_log_Passive_Sc",
                                              nullptr, nullptr, nullptr));
    block_logic[2]->SetVisAttributes(matVis2);
    // here we go again - bottome of sandwich box is
    //  -(w + sc_a + sc_p +w)/2, now add w and sc_a and half of sc_p:
    //  -w/2 -sc_a/2 - sc_p/2 -w/2 + w + sc_a + sc_p/2 = sc_a/2
    sandwich_vol[3] = new G4PVPlacement(nullptr, G4ThreeVector(0, sc_active_thickness / 2., 0),
                                        block_logic[2],
                                        "CEmc_Passive_Si_plate",
                                        sandwich, false, 0, OverlapCheck());
  }
  else
  {
    sandwich_vol[3] = nullptr;
  }
  return 0;
}

void PHG4CEmcTestBeamDetector::CalculateGeometry()
{
  // calculate the inner/outer radius of the g4tub using the test setup numbers
  // 1mm scintillator, 0.5mm tungsten, 0.25 mm gap at the back of the detector
  // the geometry with figure is explained in our wiki:
  // https://www.phenix.bnl.gov/WWW/offline/wikioff/index.php/CEmc
  // get the alpha angle via law of cos (a is the 0.25 mm gap):
  // a^2 = b^2+c^2 - 2bc*cos(alpha)
  // cos(alpha) =  (b^2+c^2 - a^2)/2bc
  // with b=c
  // cos(alpha) = 1-(a^2/2b^2)
  alpha = acos(1 - (gap * gap) / (2 * plate_x * plate_x));
  // inner radius
  inner_radius = sandwich_thickness / tan(alpha);
  outer_radius = sqrt((inner_radius + plate_x) * (inner_radius + plate_x) + sandwich_thickness * sandwich_thickness);

  // twice no_overlap to apply safety margin also for towers
  // angular coverage of 1 tower is 12 sandwiches = sandwiches_per_tower*alpha
  tower_angular_coverage = sandwiches_per_tower * alpha + 0.00005 * deg;  // safety margin added to remove 142 nm overlap
  // cemc prototype has 7 towers
  cemc_angular_coverage = num_towers * tower_angular_coverage;
  return;
}
