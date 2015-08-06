#include "PHG4OuterHcalDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeomv3.h"

#include <g4main/PHG4Utils.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <cmath>
#include <sstream>

using namespace std;

static double no_overlap = 0.00015 * cm; // added safety margin against overlaps by using same boundary between volumes

PHG4OuterHcalDetector::PHG4OuterHcalDetector( PHCompositeNode *Node, const std::string &dnam, const int lyr  ):
  PHG4Detector(Node, dnam),
  steel_rectangle_plate_x(657.2*mm),
  steel_plate_x(828.7*mm),
  steel_plate_z(3049.1*2*mm),
  n_steel_plates(320),
  scinti_tile_x(821.1*mm),
  scinti_tile_y(7*mm),
  scinti_tile_z(steel_plate_z),
  scinti_gap(8.5*mm),
  scinti_eta_coverage(1.1),
  scinti_gap_neighbor(2*mm),
  n_scinti_tiles(11),
  gap_between_tiles(2*mm),
  envelope_inner_radius(1780*mm),
  envelope_outer_radius(2603*mm),
  envelope_z(steel_plate_z*mm+no_overlap),
  tilt_angle(12*deg),
  etacutline(0.8),
  cutbox_x((steel_plate_x - steel_rectangle_plate_x*mm)*2),// twice the size we need to cut to make geo easier
  // trapezoid twice the size of what we need
  cuttrapezoid_x(cutbox_x),
  place_in_x(0 * cm),
  place_in_y(0 * cm),
  place_in_z(0 * cm),
  x_rot(0),
  y_rot(0),
  z_rot(0),
  active(0),
  absorberactive(0),
  layer(lyr),
  blackhole(0),
  steplimits(NAN),
  scintilogicnameprefix("HcalOuterScinti"),
  field_setup(NULL)
{
  double thetacutline = M_PI/2. - PHG4Utils::get_theta(etacutline);
  // okay another complication. The cut is along eta=etacutline but the box and trapezoids are calculated from the other side 
  double z_intersect = (envelope_inner_radius+(steel_plate_x-steel_rectangle_plate_x)) * tan(thetacutline);
  double z_intersect_inner = (envelope_inner_radius-(steel_plate_x-steel_rectangle_plate_x)) * tan(thetacutline);
  z_intersect = steel_plate_z/2. - z_intersect;
  z_intersect_inner = steel_plate_z/2. - z_intersect_inner;
  cutbox_z = z_intersect*2;  // twice the size we need to cut to make geo easier
  cuttrapezoid_z_short = z_intersect;
  cuttrapezoid_z_long = z_intersect+(z_intersect_inner-z_intersect);
  // inner steel plate surface:
  // scintilator gap is fixed, calculate size of the steel plate:
  // get the angle covered by one tile (360/n_steel_plates)
  // since the plate is straight we cannot just use the circumference
  // but use a straight line (surface of plate) touching the circle in its middle
  //                                  \  |
  //                         circle    \ |
  //                         (ahem)     \|
  // ------------------------------------|
  //                                    /|
  //                                   / |
  //                                  /  |
  // the length of the surface is 2*(tan(delta_phi/2)*radius)
  // then we subtract the scintilator gap from that and have the size of the steel plate
  // surface we need. Then adjust for tilt angle to get the steel_plate_yin we use
  // for creating the steel plate (it is a bit backwards, it started mistakenly with
  // using the tilted steel surface)
  G4double phi = 2*M_PI/n_steel_plates;
  G4double total_surface = 2*tan(phi/2.)*envelope_inner_radius;
  G4double steel_surface = total_surface - scinti_gap;
  // 1*mm/320. is a fudge factor, if using the calculated numbers
  // I get a 130um overlap when putting the last scintilator in
  // We have an air gap of (8.5-7)/2mm = 0.75mm + 0.13mm = 0.88mm ~ 1mm too much
  // We divide this by the number of panels to adjust the thickness of each steel plate
  G4double fudge_subtract = 1.*mm/n_steel_plates;
  steel_plate_yin = steel_surface/cos(tilt_angle/rad)-fudge_subtract;
  cutbox_y = steel_plate_yin * 10; // 10 times this size to get complete overlap even when tilted
  cuttrapezoid_y = steel_plate_yin*3; // just make it thick enough so it overlaps the steel plate in y
  // outer steel plate surface
  // scintilator gap is fixed, calculate size of the steel plate:
  // get the angle covered by one tile (360/n_steel_plates)
  // since the plate is straight we cannot just use the circumference
  // but use a straight line (surface of plate) touching the circle in its middle
  // the picture below generates a compilation error (multi line comment)
  // so this is C style comment 
  /* 
     |\
     | \  circle
     |  \ (ahem)
     ------------------------------------|
     |  /
     | /
     |/
  */
  // the length of the surface is 2*(sin(delta_phi/2)*radius)
  // then we subtract the scintilator gap from that and have the size of the steel plate
  // surface we need. Then adjust for tilt angle to get the steel_plate_yout we use
  // for creating the steel plate (it is a bit backwards, it started mistakenly with
  // using the tilted steel surface)
  total_surface = 2*sin(phi/2.)*(envelope_inner_radius+steel_plate_x);
  steel_surface = total_surface - scinti_gap;
  steel_plate_yout = steel_surface/cos(tilt_angle/rad)-fudge_subtract;

  // allocate memory for scintillator plates
  scinti_tiles_vec.assign(2*n_scinti_tiles,static_cast<G4VSolid *>(NULL));
  testbox_x[0] = steel_plate_x - steel_rectangle_plate_x*mm;
  testbox_y[0] = steel_plate_yin;
  testbox_z[0] = 1007.8*mm;
}

PHG4OuterHcalDetector::~PHG4OuterHcalDetector()
{
  if(field_setup)
    delete field_setup;
}

//_______________________________________________________________
//_______________________________________________________________
int
PHG4OuterHcalDetector::IsInOuterHcal(G4VPhysicalVolume * volume) const
{
  // G4AssemblyVolumes naming convention:
//     av_WWW_impr_XXX_YYY_ZZZ

// where:

//     WWW - assembly volume instance number
//     XXX - assembly volume imprint number
//     YYY - the name of the placed logical volume
//     ZZZ - the logical volume index inside the assembly volume
// e.g. av_1_impr_82_HcalOuterScinti_11_pv_11
// 82 the number of the scintillator mother volume
// HcalOuterScinti_11: name of scintillator slat
// 11: number of scintillator slat logical volume
  if (absorberactive)
    {
      if (steel_absorber_vec.find(volume) != steel_absorber_vec.end())
	{
	  return -1;
	}
    }
  if (volume->GetName().find(scintilogicnameprefix) != string::npos)
    {
      return 1;
    }
  return 0;
}

void
PHG4OuterHcalDetector::Construct( G4LogicalVolume* logicWorld )
{
  field_setup = new PHG4OuterHcalFieldSetup(/*G4int steelPlates*/
  n_steel_plates,
  /*G4double scintiGap*/scinti_gap, /*G4double tiltAngle*/tilt_angle);


  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4VSolid* hcal_envelope_cylinder = new G4Tubs("Hcal_envelope_solid",  envelope_inner_radius, envelope_outer_radius, envelope_z/2.,0,2*M_PI);
  G4LogicalVolume* hcal_envelope_log =  new G4LogicalVolume(hcal_envelope_cylinder, Air, G4String("Hcal_envelope"), 0, 0, 0);
  G4VisAttributes* hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Magenta());
  hcal_envelope_log->SetVisAttributes(hcalVisAtt);
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(x_rot);
  hcal_rotm.rotateY(y_rot);
  hcal_rotm.rotateZ(z_rot);
  new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(place_in_x, place_in_y, place_in_z)), hcal_envelope_log, "OuterHcal", logicWorld, 0, false, overlapcheck);
  ConstructOuterHcal(hcal_envelope_log);
  AddGeometryNode();
  return;
}

int
PHG4OuterHcalDetector::ConstructOuterHcal(G4LogicalVolume* hcalenvelope)
{
  G4VSolid *steel_plate_4 = ConstructHcalSteel(hcalenvelope);
  //   DisplayVolume(steel_plate_4 ,hcalenvelope);
  G4LogicalVolume *steel_logical = new G4LogicalVolume(steel_plate_4, G4Material::GetMaterial("G4_Fe"), "HcalOuterSteelPlate", 0, 0, 0);
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(true);
  visattchk->SetColour(G4Colour::Cyan());
  steel_logical->SetVisAttributes(visattchk);
  G4AssemblyVolume *scinti_mother_logical = ConstructHcalScintillatorAssembly(hcalenvelope);
  // visattchk = new G4VisAttributes();
  // visattchk->SetVisibility(true);
  // visattchk->SetForceSolid(false);
  // visattchk->SetColour(G4Colour::Red());
  // scinti_mother_logical->SetVisAttributes(visattchk);

  double thickness = steel_plate_yin * cos(tilt_angle / rad);
  thickness += scinti_gap;
  double deltaphi = acos((2 * envelope_inner_radius * envelope_inner_radius - thickness * thickness) / (2 * envelope_inner_radius * envelope_inner_radius));
  double phi = 0;
  ostringstream name;
  for (int i = 0; i < n_steel_plates; i++)
    {
      G4RotationMatrix *Rot = new G4RotationMatrix();
      name.str("");
      name << "HcalOuterScintiMother_" << i;
      Rot->rotateZ(phi * rad - tilt_angle);
      double xpos_scinti = (envelope_inner_radius+(scinti_tile_x)/2.) * cos(phi);
      double ypos_scinti = (envelope_inner_radius+(scinti_tile_x)/2.) * sin(phi);
      G4ThreeVector g4vec(xpos_scinti, ypos_scinti, 0);
      scinti_mother_logical->MakeImprint(hcalenvelope,g4vec,Rot,i,overlapcheck);
      Rot = new G4RotationMatrix();
      Rot->rotateZ(-phi * rad + tilt_angle);
      name.str("");
      name << "HcalOuterSteelPlate" << i;
      // start at the same position as the scintillator tiles. Naturally G4 has a different center for
      // rotating/translating the objects - the reference for the extruded solid seems to be
      // the upper left corner, the G4Box for the scinitllator has the center as reference
      // now shift this into the middle of the gap between the scintillator tiles at the inner radius
      // using the tilt angle and the angle of the slat
      double xpos = xpos_scinti;
      double ypos = ypos_scinti;
      // ypos += ((envelope_outer_radius-envelope_inner_radius)/2.)*sin(-phi/rad + tilt_angle/rad);
      // xpos -= ((envelope_outer_radius-envelope_inner_radius)/2.)*cos(-phi/rad + tilt_angle/rad);
      ypos += (scinti_tile_x/2.)*sin(-phi/rad + tilt_angle/rad);
      xpos -= (scinti_tile_x/2.)*cos(-phi/rad + tilt_angle/rad);
      // now shift the steel extruded shape down by half the width of the scintillator + gap
      // using the tilt angle and the angle of the slat
      xpos -= sin(-phi/rad + tilt_angle/rad)* (scinti_tile_y+(scinti_gap-scinti_tile_y)/2.)/2.;
      ypos -= cos(-phi/rad + tilt_angle/rad)* (scinti_tile_y+(scinti_gap-scinti_tile_y)/2.)/2.;
      steel_absorber_vec.insert(new G4PVPlacement(Rot, G4ThreeVector(xpos, ypos, 0), steel_logical, name.str().c_str(), hcalenvelope, 0, i, overlapcheck));
      phi += deltaphi;
    }

  //field after burner
  //assign the gap field strength to the air volume around HCal. Not so right for the forward corner pieces
  hcalenvelope->SetFieldManager(
      field_setup -> get_Field_Manager_Gap(),
      false);

  steel_logical->SetFieldManager(
      field_setup -> get_Field_Manager_Iron(),
      true);
//  cout <<"PHG4OuterHcalDetector::ConstructOuterHcal - register field after burner"<<endl;



  return 0;
}


G4VSolid *
PHG4OuterHcalDetector::ConstructHcalSteel(G4LogicalVolume* hcalenvelope)
{
  // ---------------------------------
  // |                               |
  // |                               |
  // |                               |
  // ---------               ---------
  //          \             /
  //           \___________/

  // strategy - use extruded solid to get the trapezoid and then
  // cut the lower parts to fit our magnet envelope
  // ---------------------------------
  // |                               |
  // |                               |
  // |                               |
  // ---------               ---------
  // |        \             /        |
  // |_________\___________/_________|
  // get corners of front piece, it is horizontal so the upper
  // corners are easy (0/0 and steel_plate_y/0)
  G4double upperleft_x = 0.; 
  G4double upperleft_y = 0.; 
  G4double upperright_x = steel_plate_x;
  G4double upperright_y = 0.;
  // we have to adjust for the tilt angle
  G4double lowerright_x = steel_plate_x+sin(tilt_angle/rad)*steel_plate_yout;
  G4double lowerright_y = -cos(tilt_angle/rad)*steel_plate_yout;
  G4double lowerleft_x = sin(tilt_angle/rad)*steel_plate_yin;
  G4double lowerleft_y = -cos(tilt_angle/rad)*steel_plate_yin;
  G4TwoVector v1( upperleft_x,upperleft_y );
  G4TwoVector v2(upperright_x,upperright_y );
  G4TwoVector v3(lowerright_x,lowerright_y );
  G4TwoVector v4(lowerleft_x,lowerleft_y );
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(v1);
  vertexes.push_back(v2);
  vertexes.push_back(v3);
  vertexes.push_back(v4);
  G4TwoVector zero(0, 0);
  G4VSolid* steel_plate =  new G4ExtrudedSolid("SteelPlate",
					       vertexes,
					       steel_plate_z  / 2.0,
					       zero, 1.0,
					       zero, 1.0);

  //  DisplayVolume(  steel_plate,hcalenvelope);
  // now cut out the lower pieces to fit this with the magnet envelope
  // ---------------------------------
  // |                               |
  // |                               |
  // |                               |
  // ---------               ---------
  //          \             /
  //           \-----------/
  // strategy - first use simple G4boxes (which need to be tilted to
  // accomodate for the 12 deg tilt_angle)
  // ---------------------------------
  // |                               |
  // |                               |
  // |                               |
  // ---------               ---------
  // |  X    |\             /|   X   |
  // |_______| \___________/ |_______|

  G4VSolid *cutbox = new G4Box("CutBox",cutbox_x/2.,cutbox_y/2.,cutbox_z/2.);
  // G4LogicalVolume* checksolid = new G4LogicalVolume(cutbox,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
  // G4VisAttributes* visattchk = new G4VisAttributes();
  // visattchk->SetVisibility(true);
  // visattchk->SetForceSolid(false);
  // visattchk->SetColour(G4Colour::Green());
  // checksolid->SetVisAttributes(visattchk);
  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateZ(-12*deg);
  double yshift = tan(tilt_angle/rad)*(steel_plate_x - steel_rectangle_plate_x*mm);
  //  new G4PVPlacement(rotm,G4ThreeVector(0,-yshift/2.,steel_plate_z/2.),checksolid,"CHKVOL",hcalenvelope, 0, false, overlapcheck);
  G4VSolid *steel_plate_1 = new G4SubtractionSolid("steel_plate_1",steel_plate,cutbox,rotm,G4ThreeVector(0,-yshift/2.,steel_plate_z/2.));
  G4VSolid *steel_plate_2 = new G4SubtractionSolid("steel_plate_2",steel_plate_1,cutbox,rotm,G4ThreeVector(0,-yshift/2.,-steel_plate_z/2.));
  //   DisplayVolume(steel_plate_2 ,hcalenvelope);

  // strategy - now cut out the leftover triangles (trapezoids, twice the size we need to make 
  // alignment easy)
  //
  // ---------               ---------|
  //          \             /         |
  //           \___________/          |
  //                      /           |
  //                     /____________|
  //
  // ---------------------------------
  // |                               |
  // |                               |
  // |                               |
  // ---------               ---------
  //          \             /
  //           \___________/

  G4VSolid *cuttrapezoid = new G4Trap("Cuttrapezoid",cuttrapezoid_y,cuttrapezoid_x,cuttrapezoid_z_long,cuttrapezoid_z_short);
  // DisplayVolume(cuttrapezoid,hcalenvelope);
  rotm = new G4RotationMatrix();
  rotm->rotateX(-90*deg);
  rotm->rotateZ(90*deg);
  G4VSolid *steel_plate_3 = new G4SubtractionSolid("steel_plate_3",steel_plate_2,cuttrapezoid,rotm,G4ThreeVector(0,0,steel_plate_z/2.-(cuttrapezoid_z_short+cuttrapezoid_z_long)/4.));
  rotm = new G4RotationMatrix();
  rotm->rotateX(90*deg);
  rotm->rotateZ(90*deg);
  G4VSolid *steel_plate_4 = new G4SubtractionSolid("steel_plate_4",steel_plate_3,cuttrapezoid,rotm,G4ThreeVector(0,0,-steel_plate_z/2.+(cuttrapezoid_z_short+cuttrapezoid_z_long)/4.));
  //   DisplayVolume(steel_plate_4 ,hcalenvelope);
  return steel_plate_4;
}

G4VSolid *
PHG4OuterHcalDetector::ConstructHcalScintillator(G4LogicalVolume* hcalenvelope)
{
  // ---------------------------------
  // |                               |
  // |                               |
  // |                               |
  // ---------               ---------
  //          \             /
  //           \___________/

  // strategy - make a G4Box for the scintillator
  // and cut out a trapezoid to fit to our magnet envelope
  // ---------------------------------
  // |                               |
  // |                               |
  // |                               |
  // ---------               ---------
  // |        \             /        |
  // |_________\___________/_________|
  G4VSolid* scinti_tile =  new G4Box("ScintiTile",scinti_tile_x/2.,scinti_tile_y/2.,scinti_tile_z/2.);

  //    DisplayVolume(  scinti_tile,hcalenvelope);
  // ---------               ---------|
  //          \             /         |
  //           \___________/          |
  //                      /           |
  //                     /____________|
  //
  // ---------------------------------
  // |                               |
  // |                               |
  // |                               |
  // ---------               ---------
  //          \             /
  //           \___________/

  G4double z_offset = 10*cm; // to make sure its size exceeds the scintillator
  G4double trap_y = scinti_tile_y*2;
  G4VSolid *cuttrapezoid = new G4Trap("Cuttrapezoid",trap_y,cuttrapezoid_x,cuttrapezoid_z_long+z_offset,cuttrapezoid_z_short+z_offset);
  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateX(-90*deg);
  rotm->rotateZ(90*deg);
  // G4LogicalVolume* checksolid = new G4LogicalVolume(cuttrapezoid,G4Material::GetMaterial("G4_POLYSTYRENE"),"DISPLAYLOGICAL", 0, 0, 0);
  // G4VisAttributes* visattchk = new G4VisAttributes();
  // visattchk->SetVisibility(true);
  // visattchk->SetForceSolid(false);
  // checksolid->SetVisAttributes(visattchk);
  // new G4PVPlacement(rotm,G4ThreeVector(-scinti_tile_x/2.,0,scinti_tile_z/2.-(cuttrapezoid_z_short+cuttrapezoid_z_long)/4.+z_offset/2.),checksolid,"DISPLAYVOL",hcalenvelope, 0, false, overlapcheck);
  //    visattchk->SetColour(G4Colour::Yellow());

  G4VSolid *scinti_tile_1 = new G4SubtractionSolid("scinti_tile_1",scinti_tile,cuttrapezoid,rotm,G4ThreeVector(-scinti_tile_x/2.,0,scinti_tile_z/2.-(cuttrapezoid_z_short+cuttrapezoid_z_long)/4.+z_offset/2.));
  rotm = new G4RotationMatrix();
  rotm->rotateX(90*deg);
  rotm->rotateZ(90*deg);
  G4VSolid *scinti_tile_2 = new G4SubtractionSolid("scinti_tile_2",scinti_tile_1,cuttrapezoid,rotm,G4ThreeVector(-scinti_tile_x/2.,0,-scinti_tile_z/2.+(cuttrapezoid_z_short+cuttrapezoid_z_long)/4.-z_offset/2.));
  return scinti_tile_2;
}

int
PHG4OuterHcalDetector::DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol,G4RotationMatrix *rotm )
{
  static int i = 0;
  G4LogicalVolume* checksolid = new G4LogicalVolume(volume,G4Material::GetMaterial("G4_POLYSTYRENE"),"DISPLAYLOGICAL", 0, 0, 0);
  G4VisAttributes* visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  switch(i)
    {
    case 0:
      visattchk->SetColour(G4Colour::Red());
      i++;
      break;
    case 1:
      visattchk->SetColour(G4Colour::Magenta());
      i++;
      break;
    case 2:
      visattchk->SetColour(G4Colour::Yellow());
      i++;
      break;
    case 3:
      visattchk->SetColour(G4Colour::Blue());
      i++;
      break;
    case 4:
      visattchk->SetColour(G4Colour::Cyan());
      i++;
      break;
    default:
      visattchk->SetColour(G4Colour::Green());
      i=0;
      break;
    }

  checksolid->SetVisAttributes(visattchk);
  new G4PVPlacement(rotm,G4ThreeVector(0,0,0),checksolid,"DISPLAYVOL",logvol, 0, false, overlapcheck);
  return 0;
}

int
PHG4OuterHcalDetector::ConstructHcalSingleScintillator(G4LogicalVolume* hcalenvelope)
{
  G4VSolid *bigtile = ConstructHcalScintillator(hcalenvelope);
  G4double delta_eta = scinti_eta_coverage/n_scinti_tiles;
  G4double offset = 10*cm;
  G4double inner_reference = envelope_inner_radius-offset;
  G4double outer_reference = envelope_outer_radius+offset;
  G4double x[4];
  G4double z[4];
  G4double eta = 0;
  G4double theta;

  for (int j=0; j<n_scinti_tiles;j++)
    {
      theta = M_PI/2 - PHG4Utils::get_theta(eta); // theta = 90 for eta=0
      x[0] = inner_reference;
      z[0] = tan(theta)*inner_reference;
      x[1] = outer_reference;
      z[1] = tan(theta)*outer_reference;
      x[2] = outer_reference;
      eta+=delta_eta;
      theta = M_PI/2 - PHG4Utils::get_theta(eta);
      z[2] = tan(theta)*outer_reference;
      x[3] = inner_reference;
      z[3] =  tan(theta)*inner_reference;
      // apply gap between scintillators
      z[0] += scinti_gap_neighbor/2.;
      z[1] += scinti_gap_neighbor/2.;
      z[2] -= scinti_gap_neighbor/2.;
      z[3] -= scinti_gap_neighbor/2.;
      vector<G4TwoVector> vertexes;
      for (int i=0; i<4; i++)
	{
	  G4TwoVector v(x[i],z[i]);
	  vertexes.push_back(v);
	}
      G4TwoVector zero(0, 0);
      G4VSolid *scinti =  new G4ExtrudedSolid("ScintillatorTile",
					      vertexes,
					      scinti_tile_y,
					      zero, 1.0,
					      zero, 1.0);
      G4RotationMatrix *rotm = new G4RotationMatrix();
      rotm->rotateX(-90*deg);
      G4VSolid *scinti_tile =  new G4IntersectionSolid("scintillator",bigtile,scinti,rotm,G4ThreeVector(-(inner_reference+outer_reference)/2., 0, 0));
      scinti_tiles_vec[j+n_scinti_tiles] = scinti_tile;
      rotm = new G4RotationMatrix();
      rotm->rotateX(90*deg);
      scinti_tile =  new G4IntersectionSolid("scintillator",bigtile,scinti,rotm,G4ThreeVector(-(inner_reference+outer_reference)/2., 0, 0));
      scinti_tiles_vec[n_scinti_tiles-j-1] =  scinti_tile;
    }
   // for (unsigned int i=0; i<scinti_tiles_vec.size(); i++)
   //   {
   //     if (scinti_tiles_vec[i])
   // 	 {
   // 	   DisplayVolume(scinti_tiles_vec[i],hcalenvelope );
   // 	 }
   //   }
  return 0;
}

G4AssemblyVolume *
PHG4OuterHcalDetector::ConstructHcalScintillatorAssembly(G4LogicalVolume* hcalenvelope)
{
  ConstructHcalSingleScintillator(hcalenvelope);
  G4AssemblyVolume *assmeblyvol = new G4AssemblyVolume();
  ostringstream name;
  G4ThreeVector g4vec;
  for (unsigned int i=0; i<scinti_tiles_vec.size(); i++)
    {
      name.str("");
      name << scintilogicnameprefix << i;
      G4UserLimits *g4userlimits = NULL;
      if (isfinite(steplimits))
	{
	  g4userlimits = new G4UserLimits(steplimits);
	}
      G4LogicalVolume *scinti_tile_logic = new G4LogicalVolume(scinti_tiles_vec[i],G4Material::GetMaterial("G4_POLYSTYRENE"),name.str().c_str(), NULL, NULL, g4userlimits);
      assmeblyvol->AddPlacedVolume(scinti_tile_logic,g4vec, NULL);

      //field after burner
      scinti_tile_logic->SetFieldManager(
          field_setup -> get_Field_Manager_Gap(),
          true);

    }
  return assmeblyvol;
}

void
PHG4OuterHcalDetector::AddGeometryNode()
{
  if (active)
    {
      ostringstream geonode;
      if (superdetector != "NONE")
	{
           geonode << "CYLINDERGEOM_" << superdetector;
        }
      else
        {
          geonode << "CYLINDERGEOM_" << detector_type << "_" << layer;
        }
      PHG4CylinderGeomContainer *geo =  findNode::getClass<PHG4CylinderGeomContainer>(topNode , geonode.str().c_str());
      if (!geo)
        {
          geo = new PHG4CylinderGeomContainer();
          PHNodeIterator iter( topNode );
          PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
          PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo, geonode.str().c_str(), "PHObject");
          runNode->addNode(newNode);
        }
      // here in the detector class we have internal units, convert to cm
      // before putting into the geom object
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeomv3(envelope_inner_radius / cm, (place_in_z - steel_plate_z / 2.) / cm, (place_in_z + steel_plate_z / 2.) / cm, (envelope_outer_radius-envelope_inner_radius) / cm, n_steel_plates,  tilt_angle/rad, 0);
      geo->AddLayerGeom(layer, mygeom);
      geo->identify();
    }
}
