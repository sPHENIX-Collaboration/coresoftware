#include "PHG4HcalTestBeamDetector.h"

#include <g4main/PHG4Utils.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Trap.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <sstream>

using namespace std;

static double no_overlap = 0.00015 * cm; // added safety margin against overlaps by using same boundary between volumes
static const double inch = 2.54 * cm;

PHG4HcalTestBeamDetector::PHG4HcalTestBeamDetector( PHCompositeNode *Node, const std::string &dnam, const int lyr  ):
  PHG4Detector(Node, dnam),
  inner_scinti(NULL),
  inner_absorber(NULL),
  outer_scinti(NULL),
  outer_absorber(NULL),
  outer_plate_x(27.1 * inch - 1 * inch), // subtract overlap
  outer_plate_z(29.53 * inch),
  outer_steel_x(outer_plate_x),
  outer_steel_y_in(1.099 * inch),
  outer_steel_y_out(1.776 * inch),
  outer_steel_z(outer_plate_z),
  inner_steel_x(12.43*inch),
  inner_steel_y_in(0.787*inch),
  inner_steel_y_out(1.115*inch),
  inner_steel_z(17.1*inch),
  place_in_x(0 * cm),
  place_in_y(0 * cm),
  place_in_z(0 * cm),
  x_rot(0),
  y_rot(0),
  z_rot(0),
  outer_tilt_angle(NAN),
  inner_tilt_angle(NAN),
  inner_radius(NAN),
  outer_radius(NAN),
  num_sandwiches(15),
  active(0),
  absorberactive(0),
  layer(lyr),
  blackhole(0)
{
  InnerGeom = new HcalGeom("INNER");
  OuterGeom = new HcalGeom("OUTER");
  outer_sc_dimension[0] =  649.8*mm;
  outer_sc_dimension[1] = 7 * mm;
  outer_sc_dimension[2] = outer_plate_z;
  inner_sc_dimension[0] =  inner_steel_x-10*mm; //4*108.6*mm; //77.4*4*mm front, back is longer it is a trapezoid 4*108.6*mm
  inner_sc_dimension[1] = 7*mm; 
  inner_sc_dimension[2] = inner_steel_z;
}

PHG4HcalTestBeamDetector::~PHG4HcalTestBeamDetector()
{
  delete InnerGeom;
  delete OuterGeom;
}

//_______________________________________________________________
//_______________________________________________________________
int
PHG4HcalTestBeamDetector::IsInHcalTestBeam(G4VPhysicalVolume * volume) const
{
  if (active && volume == inner_scinti)
    {
      return 1;
    }
  if (active && volume == outer_scinti)
    {
      return 2;
    }
  if (absorberactive && volume == inner_absorber)
    {
      return -1;
    }
  if (absorberactive && volume == outer_absorber)
    {
      return -2;
    }
  return 0;
}

void
PHG4HcalTestBeamDetector::Construct( G4LogicalVolume* logicWorld )
{
  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  CalculateGeometry();
  G4VSolid* hcal_envelope_box = new G4Box("Hcal_env",  InnerGeom->steel_plate_x/2.+OuterGeom->steel_plate_x/2.+40*cm,(OuterGeom->num_sandwiches*(OuterGeom->steel_plate_y_out+OuterGeom->scintillator_y)+OuterGeom->steel_plate_y_out)/2.+40*cm,OuterGeom->steel_plate_z/2.+1*cm);
  G4LogicalVolume* hcal_log =  new G4LogicalVolume(hcal_envelope_box, Air, G4String("Hcal_Log"), 0, 0, 0);
  G4VisAttributes* hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Magenta());
  hcal_log->SetVisAttributes(hcalVisAtt);
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(x_rot);
  hcal_rotm.rotateY(y_rot);
  hcal_rotm.rotateZ(z_rot);
  new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(place_in_x, place_in_y, place_in_z)), hcal_log, "Hcal", logicWorld, 0, false, overlapcheck);
  //   ConstructOuterHcalTower(hcal_log);
  //  ConstructOuterHcal(hcal_log);
  //  ConstructInnerHcal(hcal_log);
  //    ConstructInnerAbsorber(hcal_log);
    ConstructAbsorber(hcal_log);
    //   ConstructTest(hcal_log);
    return;
}

int
PHG4HcalTestBeamDetector::ConstructTest(G4LogicalVolume* hcalbox)
{
  G4VSolid*  inner_steel_trd = new G4Trap("Hcal_inner_steel_trd", InnerGeom->steel_plate_z, InnerGeom->steel_plate_x, InnerGeom->steel_plate_y_out, InnerGeom->steel_plate_y_in);
  //  double newx = (13.89-12.43)*inch*2;
  double tiltangle = 23.63 * deg;
  G4RotationMatrix *inner_rotm = new G4RotationMatrix;
  //inner_rotm->rotateZ(-tiltangle);
  double cutside_x = (13.89 - 12.43) * inch;
  double cutside_long =  cutside_x / cos(tiltangle / rad);
  double cutside_y = cutside_x * tan(tiltangle / rad);
  cout << "cutside_x: " << cutside_x / inch << ", cutside_longd: " << cutside_long / inch
       << ", cutside_y: " << cutside_y / inch << endl;
  G4VSolid* bigbox = new G4Trap("cutoutinner", 1.1 * InnerGeom->steel_plate_z, 1.1 * cutside_x, 1.1 * cutside_y, 0.1 * mm);
  double smalla = InnerGeom->steel_plate_y_out - cutside_y;
  double smallb = smalla * tan(tiltangle / rad);
  G4VSolid *smallbox = new G4Trap("cutsmall", 1.1 * InnerGeom->steel_plate_z, 1.1 * smalla, 1.1 * smallb, 0.1 * mm);
  //  G4SubtractionSolid *inner_steel_sub = new G4SubtractionSolid("innersub",inner_steel_trd,bigbox,NULL,G4ThreeVector(cutside_y/2.,12.43*inch+cutside_x/2.,0));
  //  G4UnionSolid *inner_steel_sub = new G4UnionSolid("innersub",inner_steel_trd,bigbox,NULL,G4ThreeVector(-cutside_y/2.,-12.43*inch/2.,0));
  G4SubtractionSolid *inner_steel_sub = new G4SubtractionSolid("innersub", inner_steel_trd, bigbox, NULL, G4ThreeVector(-cutside_y / 2., -12.43 * inch / 2., 0));
  inner_rotm->rotateZ(270 * deg);
  //  G4UnionSolid *inner_steel_sub2 = new G4UnionSolid("innersub",inner_steel_sub,smallbox,inner_rotm,G4ThreeVector(InnerGeom->steel_plate_y_in/2.,-InnerGeom->steel_plate_x/2.+(smallb+0.1*mm)/4.,0));
  G4SubtractionSolid *inner_steel_sub2 = new G4SubtractionSolid("innersub", inner_steel_sub, smallbox, inner_rotm, G4ThreeVector(InnerGeom->steel_plate_y_in / 2., -InnerGeom->steel_plate_x / 2. + (smallb + 0.1 * mm) / 4., 0));
  G4VSolid*  outer_steel_trd = new G4Trap("Hcal_outer_steel_trd", OuterGeom->steel_plate_z, OuterGeom->steel_plate_x, OuterGeom->steel_plate_y_out, OuterGeom->steel_plate_y_in);
  // double tiltangle=23.63*deg;
G4RotationMatrix *inner_rotm2 = new G4RotationMatrix;
inner_rotm2->rotateZ(tiltangle);
 double horizshift = sin(tiltangle/rad)*InnerGeom->steel_plate_x/2.-0.624*inch-InnerGeom->steel_plate_y_out/2.;
 double vertshift = cos(tiltangle/rad)*InnerGeom->steel_plate_x/2.-InnerGeom->steel_plate_x/2.;
// G4UnionSolid *outer_steel_sub = new G4UnionSolid("outerabs",outer_steel_trd,inner_steel_trd ,inner_rotm,G4ThreeVector(-(OuterGeom->steel_plate_y_out+ OuterGeom->steel_plate_y_out)/2.+(InnerGeom->steel_plate_y_out+ InnerGeom->steel_plate_y_out)/2., OuterGeom->steel_plate_x/2.+InnerGeom->steel_plate_x/2.,0));
// this shifts y so both inner and outer trapezoid line up
 double shift_y_to_zero = -(OuterGeom->steel_plate_y_out+OuterGeom->steel_plate_y_in)/4.+(InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/4.;
 double shift_x_to_zero =  OuterGeom->steel_plate_x/2.+InnerGeom->steel_plate_x/2.-(InnerGeom->steel_plate_x-12.43*inch);
 G4SubtractionSolid *outer_steel_sub3 = new G4SubtractionSolid("outerabs",outer_steel_trd,inner_steel_sub2 ,inner_rotm2,G4ThreeVector(shift_y_to_zero +horizshift ,shift_x_to_zero + vertshift ,0));
 cout << "horiz shift: " << horizshift << ", vert shift: " << vertshift << endl;
  //  G4SubtractionSolid *outer_steel_sub = new G4SubtractionSolid("outerabs",outer_steel_trd, cutout,NULL,G4ThreeVector(-(OuterGeom->steel_plate_y_out+ OuterGeom->steel_plate_y_in)/2., OuterGeom->steel_plate_x/2.,0));
  G4LogicalVolume* outer_steel_log =  new G4LogicalVolume(outer_steel_sub3, G4Material::GetMaterial("G4_Fe"), "Hcal_Outer_Steel", 0, 0, 0);
  G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetVisibility(true);
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Magenta());
  outer_steel_log->SetVisAttributes(cemcVisAtt);
  G4RotationMatrix *hcal_rotm3 = new G4RotationMatrix;
  //hcal_rotm3->rotateZ(270 * deg);
  outer_absorber = new G4PVPlacement(hcal_rotm3, G4ThreeVector( - OuterGeom->scintillator_y / 2., 0, 0), outer_steel_log, "HcalSteel", hcalbox, 0, false, overlapcheck);
  return 0;
}

int
PHG4HcalTestBeamDetector::ConstructInnerAbsorber(G4LogicalVolume* hcalbox)
{
  // first the steel plate
  G4VSolid*  inner_steel_trd = new G4Trap("Hcal_inner_steel_trd", InnerGeom->steel_plate_z, InnerGeom->steel_plate_x, InnerGeom->steel_plate_y_out, InnerGeom->steel_plate_y_in);
  //  double newx = (13.89-12.43)*inch*2;
  // now we cut out the parts to get the end right which slips into the outer steel
  double tiltangle = 23.63 * deg;
  G4RotationMatrix *inner_rotm = new G4RotationMatrix;
  //triangle to cut out
  double cutside_x = (13.89 - 12.43) * inch; // inner plate length (SP02-100-012 center drawing)
  double cutside_long =  cutside_x / cos(tiltangle / rad);
  double cutside_y = cutside_x * tan(tiltangle / rad);
  cout << "cutside_x: " << cutside_x / inch << ", cutside_longd: " << cutside_long / inch
       << ", cutside_y: " << cutside_y / inch << endl;

  G4VSolid* bigbox = new G4Trap("cutoutinner", 1.1 * InnerGeom->steel_plate_z, 1.1 * cutside_x, 1.1 * cutside_y, 0.1 * mm);
  G4SubtractionSolid *inner_steel_sub = new G4SubtractionSolid("innersub", inner_steel_trd, bigbox, NULL, G4ThreeVector(-cutside_y / 2., -12.43 * inch / 2., 0));
  double smalla = InnerGeom->steel_plate_y_out - cutside_y;
  double smallb = smalla * tan(tiltangle / rad);
  G4VSolid *smallbox = new G4Trap("cutsmall", 1.1 * InnerGeom->steel_plate_z, 1.1 * smalla, 1.1 * smallb, 0.1 * mm);
  //  G4SubtractionSolid *inner_steel_sub = new G4SubtractionSolid("innersub",inner_steel_trd,bigbox,NULL,G4ThreeVector(cutside_y/2.,12.43*inch+cutside_x/2.,0));
  //  G4UnionSolid *inner_steel_sub = new G4UnionSolid("innersub",inner_steel_trd,bigbox,NULL,G4ThreeVector(-cutside_y/2.,-12.43*inch/2.,0));
  inner_rotm->rotateZ(270 * deg);
  //  G4UnionSolid *inner_steel_sub2 = new G4UnionSolid("innersub",inner_steel_sub,smallbox,inner_rotm,G4ThreeVector(InnerGeom->steel_plate_y_in/2.,-InnerGeom->steel_plate_x/2.+(smallb+0.1*mm)/4.,0));
  G4SubtractionSolid *inner_steel_sub2 = new G4SubtractionSolid("innersub", inner_steel_sub, smallbox, inner_rotm, G4ThreeVector(InnerGeom->steel_plate_y_in / 2., -InnerGeom->steel_plate_x / 2. + (smallb + 0.1 * mm) / 4., 0));

  G4VSolid*  outer_steel_trd1 = new G4Trap("Hcal_outer_steel_trd", OuterGeom->steel_plate_z, OuterGeom->steel_plate_x, OuterGeom->steel_plate_y_out, OuterGeom->steel_plate_y_in);
  G4VSolid*  outer_steel_trd = new G4Trap("Hcal_outer_steel_trd", OuterGeom->steel_plate_z, OuterGeom->steel_plate_x, OuterGeom->steel_plate_y_out+0.7*inch, OuterGeom->steel_plate_y_in+0.7*inch);
  // double tiltangle=23.63*deg;
G4RotationMatrix *inner_rotm2 = new G4RotationMatrix;
inner_rotm2->rotateZ(tiltangle);
 double horizshift = sin(tiltangle/rad)*InnerGeom->steel_plate_x/2.-0.624*inch-((InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/2.)/2.;
 double vertshift = cos(tiltangle/rad)*InnerGeom->steel_plate_x/2.-InnerGeom->steel_plate_x/2.;
// G4UnionSolid *outer_steel_sub = new G4UnionSolid("outerabs",outer_steel_trd,inner_steel_trd ,inner_rotm,G4ThreeVector(-(OuterGeom->steel_plate_y_out+ OuterGeom->steel_plate_y_out)/2.+(InnerGeom->steel_plate_y_out+ InnerGeom->steel_plate_y_out)/2., OuterGeom->steel_plate_x/2.+InnerGeom->steel_plate_x/2.,0));
// this shifts y so both inner and outer trapezoid line up
 double shift_y_to_zero = -(OuterGeom->steel_plate_y_out+OuterGeom->steel_plate_y_in)/4.+(InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/4.;
 double shift_x_to_zero =  OuterGeom->steel_plate_x/2.+InnerGeom->steel_plate_x/2.-(InnerGeom->steel_plate_x-12.43*inch);
 G4UnionSolid *outer_steel_sub3 = new G4UnionSolid("outerabs",outer_steel_trd,inner_steel_sub2 ,inner_rotm2,G4ThreeVector(shift_y_to_zero +horizshift ,shift_x_to_zero + vertshift ,0));
 cout << "horiz shift: " << horizshift << ", vert shift: " << vertshift << endl;
  DisplayVolume(outer_steel_sub3, hcalbox);
  DisplayVolume(outer_steel_trd1, hcalbox);

  G4LogicalVolume* checksolid = new G4LogicalVolume(inner_steel_sub2,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
  G4VisAttributes* visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
      visattchk->SetColour(G4Colour::Blue());
  checksolid->SetVisAttributes(visattchk);
  G4RotationMatrix *inner_rotmtst = new G4RotationMatrix(*inner_rotm2);
  new G4PVPlacement(inner_rotmtst,G4ThreeVector(shift_y_to_zero +horizshift ,shift_x_to_zero + vertshift ,0),checksolid,"CHKVOL",hcalbox, 0, false, overlapcheck);
  //  DisplayVolume(inner_steel_sub2, hcalbox);
  return 0;
  //  G4SubtractionSolid *outer_steel_sub = new G4SubtractionSolid("outerabs",outer_steel_trd, cutout,NULL,G4ThreeVector(-(OuterGeom->steel_plate_y_out+ OuterGeom->steel_plate_y_in)/2., OuterGeom->steel_plate_x/2.,0));
  G4LogicalVolume* outer_steel_log =  new G4LogicalVolume(outer_steel_sub3, G4Material::GetMaterial("G4_Fe"), "Hcal_Outer_Steel", 0, 0, 0);
  G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetVisibility(true);
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Magenta());
  outer_steel_log->SetVisAttributes(cemcVisAtt);
  G4VSolid *OuterScinti = new G4Box ("outerscinti", 7*mm/2., OuterGeom->steel_plate_x/2.-1*inch, OuterGeom->steel_plate_z/2.);
  G4LogicalVolume* outer_scinti = new G4LogicalVolume(OuterScinti,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
  G4VisAttributes* cemcVisAtt1 = new G4VisAttributes();
  cemcVisAtt1->SetVisibility(true);
  cemcVisAtt1->SetForceSolid(false);
  cemcVisAtt1->SetColour(G4Colour::Blue());
outer_scinti->SetVisAttributes(cemcVisAtt1);
  ostringstream inner_hcal_name;
  double yup = -((InnerGeom->num_sandwiches-1)*8*mm+InnerGeom->steel_plate_y_out)/2.;
  double xup = 0.;
  //  double xoffset = -OuterGeom->steel_plate_x/2. - InnerGeom->steel_plate_x/2.;
  //return 0;
  //for (int i=0; i<InnerGeom->num_sandwiches; i++)
    for (int i=0; i<3; i++)
    {
      inner_hcal_name.str("");
      inner_hcal_name << "InnerHcal_" << i;

  G4RotationMatrix *hcal_rotm3 = new G4RotationMatrix;
  hcal_rotm3->rotateZ(270 * deg-i*InnerGeom->tilt_angle);
  //  outer_absorber = new G4PVPlacement(hcal_rotm3, G4ThreeVector( - OuterGeom->scintillator_y / 2., 0, 0), outer_steel_log, "HcalSteel", hcalbox, 0, false, overlapcheck);
  outer_absorber = new G4PVPlacement(hcal_rotm3, G4ThreeVector(xup, yup, 0), outer_steel_log, "HcalSteel", hcalbox, 0, false, overlapcheck);
  G4RotationMatrix *hcal_rotm4 = new G4RotationMatrix;
  hcal_rotm4->rotateZ(270 * deg-i*InnerGeom->tilt_angle);
  double yscinit_up = yup - (OuterGeom->steel_plate_y_in)+OuterGeom->scinti_gap;//+OuterGeom->steel_plate_y_out)/2.;
new G4PVPlacement(hcal_rotm4, G4ThreeVector(xup+1*inch, yscinit_up, 0),outer_scinti , "HcalScinti", hcalbox, 0, false, overlapcheck);
      yup += ((InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/2.+InnerGeom->scintillator_y_vertical_projection)*cos(i*InnerGeom->tilt_angle) + 2*8*mm;
      xup -=  ((InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/2.+InnerGeom->scintillator_y_vertical_projection)*sin(i*InnerGeom->tilt_angle);
    }
    return 0;
}

int
PHG4HcalTestBeamDetector::ConstructAbsorber(G4LogicalVolume* hcalbox)
{
  // first the inner steel plate
  G4VSolid*  inner_steel_trd = new G4Trap("Hcal_inner_steel_trd", InnerGeom->steel_plate_z, InnerGeom->steel_plate_x, InnerGeom->steel_plate_y_out, InnerGeom->steel_plate_y_in);
  // now we cut out the parts to get the end right which slips into the outer steel
  // first cut - remove triangle at the back bottom
  //___________________
  //                   |
  //                   |
  //                  /|
  //     23.63 deg   / | <-- G4 has no triangles so use trapezoid, upper "pointy" size 0.1mm 
  // _______________/__|
  double tiltangle = 23.63 * deg;
  //triangle to cut out
  double cutside_x = (13.89 - 12.43) * inch; // inner plate length (SP02-100-012 center drawing)
  double cutside_long =  cutside_x / cos(tiltangle / rad);
  double cutside_y = cutside_x * tan(tiltangle / rad);
  cout << "cutside_x: " << cutside_x / inch << ", cutside_longd: " << cutside_long / inch
       << ", cutside_y: " << cutside_y / inch << endl;
  // we need t osubtract a larger triangle to avoid edge effects, the triangle is multiplied by
  // a marginfactor but currently this only works for  marginfactor = 2. 
  // need to understand the general shift when doing the G4SubtractionSolid
  double marginfactor = 2.;
  G4VSolid* bigtriangle = new G4Trap("cutoutinner", 1.1 * InnerGeom->steel_plate_z, marginfactor * cutside_x,  marginfactor * cutside_y, 0.0001 * mm);
  // shifts to be applied for original triangle we want to subtract
 double yshift = cutside_y;
 double xshift = InnerGeom->steel_plate_x - cutside_x;
 // the marginfactor needs to be reworked if this should work for the general case, this only works
 // for marginfactor = 2.
 G4VSolid* inner_steel_sub = new G4SubtractionSolid("innersub", inner_steel_trd, bigtriangle, NULL, G4ThreeVector(-(yshift+(cutside_y/marginfactor)/2.) / 2., -(xshift+(cutside_x/marginfactor)/2.) / 2., 0));

  // second cut - remove smaller triangle at the top
  //___________________
  //                 \ | <-- again G4 has no triangles so use trapezoid, upper "pointy" size 0.1mm 
  //                  \|
  //                  /
  //     23.63 deg   /  
  // _______________/
  double smalla = InnerGeom->steel_plate_y_out - cutside_y;
  double smallb = smalla * tan(tiltangle / rad);
  G4VSolid *smalltriangle = new G4Trap("cutsmall", 1.1 * InnerGeom->steel_plate_z, marginfactor * smalla, marginfactor * smallb, 0.0001 * mm);
  // need to rotate small box by 270 deg before cutting it out
  G4RotationMatrix *smalltriangle_rotm = new G4RotationMatrix;
  smalltriangle_rotm->rotateZ(270 * deg);
  // same as above, the marginfactor needs to be reworked if this should work for the 
  // general case, this only works for marginfactor = 2. It probably has something
  // to do with these triangles being trapezoids with a center of a+b/2 where b=0.0001mm

    G4VSolid *inner_steel = new G4SubtractionSolid("innersub", inner_steel_sub, smalltriangle, smalltriangle_rotm, G4ThreeVector(((-(InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/2.)/2.+cutside_y+smalla/4.), (-InnerGeom->steel_plate_x / 2. - (smallb/marginfactor) / 2.), 0));
  // done with inner steel, the cutout to fit the outer steel is missing but we use a 
  // G4UnionSolid to merge the outer and inner steel so this doesn't matter and the 
  // gaps between the outer and inner steel plates are really marginal. 

  G4VSolid*  outer_steel_trd = new G4Trap("Hcal_outer_steel_trd", OuterGeom->steel_plate_z, OuterGeom->steel_plate_x, OuterGeom->steel_plate_y_out, OuterGeom->steel_plate_y_in);
  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateZ(-(90*deg-tiltangle));

  //  DisplayVolume(inner_steel, hcalbox, rotm);
  //  DisplayVolume(outer_steel_trd, hcalbox);
G4LogicalVolume* checkinnersteel = new G4LogicalVolume(inner_steel,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
  G4VisAttributes* visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
      visattchk->SetColour(G4Colour::Green());
  visattchk->SetForceSolid(false);
  checkinnersteel->SetVisAttributes(visattchk);
 // double horizshift = sin(tiltangle/rad)*InnerGeom->steel_plate_x/2.;
 // double vertshift = cos(tiltangle/rad)*InnerGeom->steel_plate_x/2.-InnerGeom->steel_plate_x/2.;
 // double shift_y_to_zero = -(OuterGeom->steel_plate_y_out+OuterGeom->steel_plate_y_in)/4.+(InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/4.;
 // shift_x_to_zero =  OuterGeom->steel_plate_x/2.+InnerGeom->steel_plate_x/2.;

//  new G4PVPlacement(rotm,G4ThreeVector(shift_y_to_zero+horizshift,shift_x_to_zero+vertshift,0),checkinnersteel,"CHKVOL",hcalbox, 0, false, overlapcheck);
  new G4PVPlacement(rotm,G4ThreeVector(0,0,0),checkinnersteel,"CHKVOL",hcalbox, 0, false, overlapcheck);
  double inner_steel_x_shift = (-OuterGeom->steel_plate_x/2.) - sin(tiltangle/rad)*InnerGeom->steel_plate_x/2.;
    //(-OuterGeom->steel_plate_x/2.)-(InnerGeom->steel_plate_x/2.);
  double inner_steel_y_shift = cos(tiltangle/rad)*InnerGeom->steel_plate_x/2./2.;
  //-(OuterGeom->steel_plate_y_in+OuterGeom->steel_plate_y_out)/2./2.;
  new G4PVPlacement(rotm,G4ThreeVector(inner_steel_x_shift,inner_steel_y_shift,0),checkinnersteel,"CHKVOL",hcalbox, 0, false, overlapcheck);
// inner_steel_y_shift += 
//   new G4PVPlacement(rotm,G4ThreeVector(inner_steel_x_shift,inner_steel_y_shift,0),checkinnersteel,"CHKVOL",hcalbox, 0, false, overlapcheck);
  //  new G4PVPlacement(rotm,G4ThreeVector(-OuterGeom->steel_plate_x/2.+1*inch/2.,-(OuterGeom->steel_plate_y_in+OuterGeom->steel_plate_y_out)/2./2.-0.624*inch/2.,0),checkinnersteel,"CHKVOL",hcalbox, 0, false, overlapcheck);

G4LogicalVolume* checksolid = new G4LogicalVolume(outer_steel_trd,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
  visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
      visattchk->SetColour(G4Colour::Magenta());
  visattchk->SetForceSolid(false);
  checksolid->SetVisAttributes(visattchk);
 G4RotationMatrix *hcal_rotm = new G4RotationMatrix;
hcal_rotm->rotateZ(270 * deg);
//  new G4PVPlacement(hcal_rotm,G4ThreeVector(5.625*inch,-2.946*inch+0.624*inch,0),checksolid,"CHKVOL",hcalbox, 0, false, overlapcheck);
 new G4PVPlacement(hcal_rotm,G4ThreeVector(0,0,0),checksolid,"CHKVOL",hcalbox, 0, false, overlapcheck);
 G4VSolid *box = new G4Box("tstbox",0.654*inch/2.,8*mm/2.,OuterGeom->scintillator_z / 2.);
checksolid = new G4LogicalVolume(box,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
// G4LogicalVolume *checksolid2 = new G4LogicalVolume(smalltriangle,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
// G4LogicalVolume *checksolid3 = new G4LogicalVolume(bigtriangle,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
 visattchk = new G4VisAttributes();
   visattchk->SetVisibility(true);
       visattchk->SetColour(G4Colour::Yellow());
   visattchk->SetForceSolid(false); 
  checksolid->SetVisAttributes(visattchk);
  new G4PVPlacement(NULL,G4ThreeVector(-OuterGeom->steel_plate_x/2.,-(OuterGeom->steel_plate_y_in+OuterGeom->steel_plate_y_out)/2,0),checksolid,"CHKVOL",hcalbox, 0, false, overlapcheck);
// visattchk = new G4VisAttributes();
//   visattchk->SetVisibility(true);
//       visattchk->SetColour(G4Colour::Magenta());
//   visattchk->SetForceSolid(false);
//  checksolid2->SetVisAttributes(visattchk);
// visattchk = new G4VisAttributes();
//   visattchk->SetVisibility(true);
//       visattchk->SetColour(G4Colour::Cyan());
//   visattchk->SetForceSolid(false);
//  checksolid3->SetVisAttributes(visattchk);
// G4LogicalVolume *isteel =  new G4LogicalVolume(inner_steel_trd,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
// visattchk = new G4VisAttributes();
//   visattchk->SetVisibility(true);
//       visattchk->SetColour(G4Colour::Green());
//   visattchk->SetForceSolid(false);
// isteel->SetVisAttributes(visattchk);
//  new G4PVPlacement(NULL,G4ThreeVector(0,0,0),isteel,"CHKVOL",hcalbox, 0, false, overlapcheck);
  //  new G4PVPlacement(rotm,G4ThreeVector(-cutside_y / 2., -12.43 * inch / 2., 0),checksolid,"CHKVOL",isteel, 0, false, overlapcheck);

  //  DisplayVolume(box1, hcalbox);
  //  DisplayVolume(outer_steel_trd1, hcalbox);
  return 0;
  outer_steel_trd = new G4Trap("Hcal_outer_steel_trd", OuterGeom->steel_plate_z, OuterGeom->steel_plate_x, OuterGeom->steel_plate_y_out+0.7*inch, OuterGeom->steel_plate_y_in+0.7*inch);
  // double tiltangle=23.63*deg;
G4RotationMatrix *inner_rotm2 = new G4RotationMatrix;
inner_rotm2->rotateZ(tiltangle);
 double horizshift = sin(tiltangle/rad)*InnerGeom->steel_plate_x/2.-0.624*inch-((InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/2.)/2.;
 double vertshift = cos(tiltangle/rad)*InnerGeom->steel_plate_x/2.-InnerGeom->steel_plate_x/2.;
// G4UnionSolid *outer_steel_sub = new G4UnionSolid("outerabs",outer_steel_trd,inner_steel_trd ,inner_rotm,G4ThreeVector(-(OuterGeom->steel_plate_y_out+ OuterGeom->steel_plate_y_out)/2.+(InnerGeom->steel_plate_y_out+ InnerGeom->steel_plate_y_out)/2., OuterGeom->steel_plate_x/2.+InnerGeom->steel_plate_x/2.,0));
// this shifts y so both inner and outer trapezoid line up
 double shift_y_to_zero = -(OuterGeom->steel_plate_y_out+OuterGeom->steel_plate_y_in)/4.+(InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/4.;
 double shift_x_to_zero =  OuterGeom->steel_plate_x/2.+InnerGeom->steel_plate_x/2.-(InnerGeom->steel_plate_x-12.43*inch);
 G4UnionSolid *outer_steel_sub3 = new G4UnionSolid("outerabs",outer_steel_trd,inner_steel ,inner_rotm2,G4ThreeVector(shift_y_to_zero +horizshift ,shift_x_to_zero + vertshift ,0));
 cout << "horiz shift: " << horizshift << ", vert shift: " << vertshift << endl;
  DisplayVolume(outer_steel_sub3, hcalbox);

  // G4LogicalVolume* checksolid = new G4LogicalVolume(inner_steel,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
  // G4VisAttributes* visattchk = new G4VisAttributes();
  // visattchk->SetVisibility(true);
  // visattchk->SetForceSolid(false);
  //     visattchk->SetColour(G4Colour::Blue());
  // checksolid->SetVisAttributes(visattchk);
  // G4RotationMatrix *inner_rotmtst = new G4RotationMatrix(*inner_rotm2);
  // new G4PVPlacement(inner_rotmtst,G4ThreeVector(shift_y_to_zero +horizshift ,shift_x_to_zero + vertshift ,0),checksolid,"CHKVOL",hcalbox, 0, false, overlapcheck);
  //  DisplayVolume(inner_steel_sub2, hcalbox);
  return 0;
  //  G4SubtractionSolid *outer_steel_sub = new G4SubtractionSolid("outerabs",outer_steel_trd, cutout,NULL,G4ThreeVector(-(OuterGeom->steel_plate_y_out+ OuterGeom->steel_plate_y_in)/2., OuterGeom->steel_plate_x/2.,0));
  G4LogicalVolume* outer_steel_log =  new G4LogicalVolume(outer_steel_sub3, G4Material::GetMaterial("G4_Fe"), "Hcal_Outer_Steel", 0, 0, 0);
  G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetVisibility(true);
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Magenta());
  outer_steel_log->SetVisAttributes(cemcVisAtt);
  G4VSolid *OuterScinti = new G4Box ("outerscinti", 7*mm/2., OuterGeom->steel_plate_x/2.-1*inch, OuterGeom->steel_plate_z/2.);
  G4LogicalVolume* outer_scinti = new G4LogicalVolume(OuterScinti,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
  G4VisAttributes* cemcVisAtt1 = new G4VisAttributes();
  cemcVisAtt1->SetVisibility(true);
  cemcVisAtt1->SetForceSolid(false);
  cemcVisAtt1->SetColour(G4Colour::Blue());
outer_scinti->SetVisAttributes(cemcVisAtt1);
  ostringstream inner_hcal_name;
  double yup = -((InnerGeom->num_sandwiches-1)*8*mm+InnerGeom->steel_plate_y_out)/2.;
  double xup = 0.;
  //  double xoffset = -OuterGeom->steel_plate_x/2. - InnerGeom->steel_plate_x/2.;
  //return 0;
  //for (int i=0; i<InnerGeom->num_sandwiches; i++)
    for (int i=0; i<3; i++)
    {
      inner_hcal_name.str("");
      inner_hcal_name << "InnerHcal_" << i;

  G4RotationMatrix *hcal_rotm3 = new G4RotationMatrix;
  hcal_rotm3->rotateZ(270 * deg-i*InnerGeom->tilt_angle);
  //  outer_absorber = new G4PVPlacement(hcal_rotm3, G4ThreeVector( - OuterGeom->scintillator_y / 2., 0, 0), outer_steel_log, "HcalSteel", hcalbox, 0, false, overlapcheck);
  outer_absorber = new G4PVPlacement(hcal_rotm3, G4ThreeVector(xup, yup, 0), outer_steel_log, "HcalSteel", hcalbox, 0, false, overlapcheck);
  G4RotationMatrix *hcal_rotm4 = new G4RotationMatrix;
  hcal_rotm4->rotateZ(270 * deg-i*InnerGeom->tilt_angle);
  double yscinit_up = yup - (OuterGeom->steel_plate_y_in)+OuterGeom->scinti_gap;//+OuterGeom->steel_plate_y_out)/2.;
new G4PVPlacement(hcal_rotm4, G4ThreeVector(xup+1*inch, yscinit_up, 0),outer_scinti , "HcalScinti", hcalbox, 0, false, overlapcheck);
      yup += ((InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/2.+InnerGeom->scintillator_y_vertical_projection)*cos(i*InnerGeom->tilt_angle) + 2*8*mm;
      xup -=  ((InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/2.+InnerGeom->scintillator_y_vertical_projection)*sin(i*InnerGeom->tilt_angle);
    }
    return 0;
}


int
PHG4HcalTestBeamDetector::ConstructOuterHcalTower(G4LogicalVolume* hcalbox)
{
  G4VSolid*  outer_steel_trd = new G4Trap("Hcal_outer_steel_trd", OuterGeom->steel_plate_z, OuterGeom->steel_plate_x, OuterGeom->steel_plate_y_out, OuterGeom->steel_plate_y_in);
  G4VSolid*  inner_steel_trd = new G4Trap("Hcal_inner_steel_trd", inner_steel_z, inner_steel_x, inner_steel_y_out, inner_steel_y_in);
 double tiltangle=23.63*deg;
G4RotationMatrix *inner_rotm = new G4RotationMatrix;
inner_rotm->rotateZ(tiltangle);
 double horizshift = sin(tiltangle/rad)*InnerGeom->steel_plate_x/2.;
 double vertshift = cos(tiltangle/rad)*InnerGeom->steel_plate_x/2.-InnerGeom->steel_plate_x/2.;
// G4UnionSolid *outer_steel_sub = new G4UnionSolid("outerabs",outer_steel_trd,inner_steel_trd ,inner_rotm,G4ThreeVector(-(OuterGeom->steel_plate_y_out+ OuterGeom->steel_plate_y_out)/2.+(InnerGeom->steel_plate_y_out+ InnerGeom->steel_plate_y_out)/2., OuterGeom->steel_plate_x/2.+InnerGeom->steel_plate_x/2.,0));
// this shifts y so both inner and outer trapezoid line up
 double shift_y_to_zero = -(OuterGeom->steel_plate_y_out+OuterGeom->steel_plate_y_in)/4.+(InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/4.;
 double shift_x_to_zero =  OuterGeom->steel_plate_x/2.+InnerGeom->steel_plate_x/2.;
 G4UnionSolid *outer_steel_sub = new G4UnionSolid("outerabs",outer_steel_trd,inner_steel_trd ,inner_rotm,G4ThreeVector(shift_y_to_zero+horizshift,shift_x_to_zero+vertshift ,0));
  //  G4SubtractionSolid *outer_steel_sub = new G4SubtractionSolid("outerabs",outer_steel_trd, cutout,NULL,G4ThreeVector(-(OuterGeom->steel_plate_y_out+ OuterGeom->steel_plate_y_in)/2., OuterGeom->steel_plate_x/2.,0));
   G4LogicalVolume* outer_steel_log =  new G4LogicalVolume(outer_steel_sub, G4Material::GetMaterial("G4_Fe"), "Hcal_Outer_Steel", 0, 0, 0);
   G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetVisibility(true);
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Magenta());
  outer_steel_log->SetVisAttributes(cemcVisAtt);
 G4RotationMatrix *hcal_rotm = new G4RotationMatrix;
			       // hcal_rotm->rotateZ(270 * deg);
outer_absorber = new G4PVPlacement(hcal_rotm, G4ThreeVector( - OuterGeom->scintillator_y / 2., 0, 0), outer_steel_log, "HcalSteel", hcalbox, 0, false, overlapcheck);
  return 0;
}

// here we build up the outer hcal from the sandwiches (steel + scintillator)
int
PHG4HcalTestBeamDetector::ConstructOuterHcal( G4LogicalVolume* logicWorld )
{
  // fudge factors to remove overlaps
  double addenvelope = 0.003*mm;
  double addverticaldistance = 0.008*mm;
  double addtoyup = 0.003*mm;

  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  double y_out = OuterGeom->steel_plate_y_out +  OuterGeom->scintillator_y_vertical_projection+addenvelope;
  double y_in = OuterGeom->steel_plate_y_in+ OuterGeom->scintillator_y_vertical_projection+addenvelope;
  G4VSolid*  outer_sandwich = new G4Trap("Hcal_Outer_Sandwich_Trap", OuterGeom->steel_plate_z,  OuterGeom->steel_plate_x,  y_out, y_in);
  G4LogicalVolume* outer_sandwich_log =  new G4LogicalVolume(outer_sandwich, Air, "Hcal_Outer_Sandwich", 0, 0, 0);
  G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Blue());
  outer_sandwich_log->SetVisAttributes(cemcVisAtt);
  ostringstream outer_hcal_name;
  double yup = -((OuterGeom->num_sandwiches-1)*y_out+OuterGeom->steel_plate_y_out)/2.;
  double xup = 0.;
  for (int i=0; i<OuterGeom->num_sandwiches; i++)
    //  for (int i=0; i<1; i++)
    {
      outer_hcal_name.str("");
      outer_hcal_name << "OuterHcal_" << i;
      G4RotationMatrix *hcal_rotm = new G4RotationMatrix;

      //      hcal_rotm->rotateZ((270-23.63/2.) * deg-i*outer_tilt_angle*rad);
      hcal_rotm->rotateZ(270 * deg - i*OuterGeom->tilt_angle);
      new G4PVPlacement(hcal_rotm, G4ThreeVector(xup, yup, 0), outer_sandwich_log,outer_hcal_name.str().c_str() , logicWorld, 0, i, overlapcheck);
      yup += ((OuterGeom->steel_plate_y_out+OuterGeom->steel_plate_y_in)/2.+OuterGeom->scintillator_y_vertical_projection+addverticaldistance)*cos(i*OuterGeom->tilt_angle)+addtoyup;
      xup -=  ((OuterGeom->steel_plate_y_out+OuterGeom->steel_plate_y_in)/2.+OuterGeom->scintillator_y_vertical_projection+addverticaldistance)*sin(i*OuterGeom->tilt_angle);
    }
  ConstructOuterSandwichVolume(outer_sandwich_log);
  // add the top absorber
  G4VSolid*  outer_steel_trd = new G4Trap("Hcal_outer_steel_trd", OuterGeom->steel_plate_z, OuterGeom->steel_plate_x, OuterGeom->steel_plate_y_out, OuterGeom->steel_plate_y_in);
  G4LogicalVolume* outer_steel_log =  new G4LogicalVolume(outer_steel_trd, G4Material::GetMaterial("G4_Fe"), "Hcal_Outer_Steel", 0, 0, 0);
  G4VisAttributes* cemcVisAtt2 = new G4VisAttributes();
  cemcVisAtt2->SetVisibility(true);
  cemcVisAtt2->SetForceSolid(false);
  cemcVisAtt2->SetColour(G4Colour::Magenta());
  outer_steel_log->SetVisAttributes(cemcVisAtt2);
  outer_hcal_name.str("");
  outer_hcal_name << "OuterHcal_" <<  OuterGeom->num_sandwiches;
  G4RotationMatrix *hcal_rotm = new G4RotationMatrix;

  //      hcal_rotm->rotateZ((270-23.63/2.) * deg-i*outer_tilt_angle*rad);
  hcal_rotm->rotateZ(270 * deg - OuterGeom->num_sandwiches * OuterGeom->tilt_angle);
  new G4PVPlacement(hcal_rotm, G4ThreeVector(xup, yup, 0), outer_steel_log, outer_hcal_name.str().c_str() , logicWorld, 0, OuterGeom->num_sandwiches, overlapcheck);

  return 0;
}

// here we build up the outer hcal from the sandwiches (steel + scintillator)
int
PHG4HcalTestBeamDetector::ConstructInnerHcal( G4LogicalVolume* logicWorld )
{
  // fudge factors to remove overlaps
  double addenvelope = 0.00*mm;
  double addverticaldistance = 0.00*mm;
  double addtoyup = 0.00*mm;

  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  double y_out = InnerGeom->steel_plate_y_out +  InnerGeom->scintillator_y_vertical_projection+addenvelope;
  double y_in = InnerGeom->steel_plate_y_in+ InnerGeom->scintillator_y_vertical_projection+addenvelope;
  G4VSolid*  inner_sandwich = new G4Trap("Hcal_Inner_Sandwich_Trap", InnerGeom->steel_plate_z,  InnerGeom->steel_plate_x,  y_out, y_in);
  G4LogicalVolume* inner_sandwich_log =  new G4LogicalVolume(inner_sandwich, Air, "Hcal_Inner_Sandwich", 0, 0, 0);
  G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Blue());
  inner_sandwich_log->SetVisAttributes(cemcVisAtt);
  ostringstream inner_hcal_name;
  double yup = -((InnerGeom->num_sandwiches-1)*y_out+InnerGeom->steel_plate_y_out)/2.;
  double xup = 0.;
  double xoffset = -OuterGeom->steel_plate_x/2. - InnerGeom->steel_plate_x/2.;
  //  for (int i=0; i<InnerGeom->num_sandwiches; i++)
   for (int i=0; i<2; i++)
    {
      inner_hcal_name.str("");
      inner_hcal_name << "InnerHcal_" << i;
      G4RotationMatrix *hcal_rotm = new G4RotationMatrix;

      hcal_rotm->rotateZ((270+23.63) * deg-i*InnerGeom->tilt_angle);
      //      hcal_rotm->rotateZ(270 * deg - i*InnerGeom->tilt_angle);
      new G4PVPlacement(hcal_rotm, G4ThreeVector(xup+xoffset, yup, 0), inner_sandwich_log,inner_hcal_name.str().c_str() , logicWorld, 0, i, overlapcheck);
      yup += ((InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/2.+InnerGeom->scintillator_y_vertical_projection+addverticaldistance)*cos(i*InnerGeom->tilt_angle)+addtoyup;
      xup -=  ((InnerGeom->steel_plate_y_out+InnerGeom->steel_plate_y_in)/2.+InnerGeom->scintillator_y_vertical_projection+addverticaldistance)*sin(i*InnerGeom->tilt_angle);
    }
  return 0;
  ConstructInnerSandwichVolume(inner_sandwich_log);
  // add the top absorber
  G4VSolid*  inner_steel_trd = new G4Trap("Hcal_inner_steel_trd", InnerGeom->steel_plate_z, InnerGeom->steel_plate_x, InnerGeom->steel_plate_y_out, InnerGeom->steel_plate_y_in);
  G4LogicalVolume* inner_steel_log =  new G4LogicalVolume(inner_steel_trd, G4Material::GetMaterial("G4_Fe"), "Hcal_Inner_Steel", 0, 0, 0);
  G4VisAttributes* cemcVisAtt2 = new G4VisAttributes();
  cemcVisAtt2->SetVisibility(true);
  cemcVisAtt2->SetForceSolid(false);
  cemcVisAtt2->SetColour(G4Colour::Magenta());
  inner_steel_log->SetVisAttributes(cemcVisAtt2);
  inner_hcal_name.str("");
  inner_hcal_name << "InnerHcal_" <<  InnerGeom->num_sandwiches;
  G4RotationMatrix *hcal_rotm = new G4RotationMatrix;

  hcal_rotm->rotateZ((270+23.63) * deg-InnerGeom->num_sandwiches*InnerGeom->tilt_angle);
  //  hcal_rotm->rotateZ(270 * deg - InnerGeom->num_sandwiches * InnerGeom->tilt_angle);
  new G4PVPlacement(hcal_rotm, G4ThreeVector(xup+xoffset, yup, 0), inner_steel_log, inner_hcal_name.str().c_str() , logicWorld, 0, InnerGeom->num_sandwiches, overlapcheck);

  return 0;
}


int
PHG4HcalTestBeamDetector::ConstructInnerHcalTower(G4LogicalVolume* logicWorld)
{
//   return 0;
// }

// int
// PHG4HcalTestBeamDetector::ConstructInnerHcal( G4LogicalVolume* logicWorld )
// {
  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4VSolid*  inner_sandwich = new G4Trap("Hcal_Inner_Sandwich_Trap", inner_steel_z, inner_steel_x, inner_steel_y_out+ inner_sc_dimension[1]+4*no_overlap, inner_steel_y_in+ inner_sc_dimension[1]+4*no_overlap);
  G4LogicalVolume* inner_sandwich_log =  new G4LogicalVolume(inner_sandwich, Air, "Hcal_Inner_Sandwich", 0, 0, 0);
  G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Blue());
inner_sandwich_log->SetVisAttributes(cemcVisAtt);
  ostringstream inner_hcal_name;
  //  double yup = -((num_sandwiches-1)*(inner_steel_y_out+inner_sc_dimension[1])+inner_steel_y_out)/2.;
  double yup = -((num_sandwiches-1)*(outer_steel_y_out+outer_sc_dimension[1])+outer_steel_y_out)/2.;
  double xup = -outer_steel_x/2.-inner_steel_x/2. - 1*inch;
  for (int i=0; i<15; i++)
    {
      inner_hcal_name.str("");
      inner_hcal_name << "InnerHcal_" << i;
      G4RotationMatrix *hcal_rotm = new G4RotationMatrix;

      hcal_rotm->rotateZ((270+23.63/2.) * deg-i*inner_tilt_angle*rad);

      new G4PVPlacement(hcal_rotm, G4ThreeVector(xup, yup, 0), inner_sandwich_log,inner_hcal_name.str().c_str() , logicWorld, 0, i, overlapcheck);
      yup += ((inner_steel_y_out+inner_steel_y_in)/2.+inner_sc_dimension[1])*cos(i*inner_tilt_angle)+0.6*mm;
      xup -=  ((inner_steel_y_out+inner_steel_y_in)/2.+inner_sc_dimension[1])*sin(i*inner_tilt_angle)+0.4*mm;
    }
  ConstructInnerSandwichVolume(inner_sandwich_log);

  return 0;
}

// here we put a single iron + scintillator sandwich together
int
PHG4HcalTestBeamDetector::ConstructOuterSandwichVolume(G4LogicalVolume* sandwich_log)
{
  vector<G4LogicalVolume*> block_logic;
  G4VSolid*  outer_steel_trd = new G4Trap("Hcal_outer_steel_trd", OuterGeom->steel_plate_z, OuterGeom->steel_plate_x, OuterGeom->steel_plate_y_out, OuterGeom->steel_plate_y_in);
  cout << "steel: z: " << OuterGeom->steel_plate_z << ", x: " << OuterGeom->steel_plate_x
       << ", yout: " <<  OuterGeom->steel_plate_y_out << ", yin: " << OuterGeom->steel_plate_y_in
       << endl;
  G4LogicalVolume* outer_steel_log =  new G4LogicalVolume(outer_steel_trd, G4Material::GetMaterial("G4_Fe"), "Hcal_Outer_Steel", 0, 0, 0);
  G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetVisibility(true);
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Magenta());
  outer_steel_log->SetVisAttributes(cemcVisAtt);
  G4VSolid *outer_scintillator = new G4Box("Outer_Scinti_Box", OuterGeom->scintillator_y/ 2., OuterGeom->scintillator_x / 2., OuterGeom->scintillator_z / 2.);
  G4LogicalVolume* outer_scinti_log = new G4LogicalVolume(outer_scintillator, G4Material::GetMaterial("G4_POLYSTYRENE"), "Hcal_Outer_Sc", 0, 0, 0);
  G4VisAttributes* cemcVisAtt1 = new G4VisAttributes();
  cemcVisAtt1->SetVisibility(true);
  cemcVisAtt1->SetForceSolid(false);
  cemcVisAtt1->SetColour(G4Colour::Green());
  outer_scinti_log->SetVisAttributes(cemcVisAtt1);
  outer_absorber = new G4PVPlacement(NULL, G4ThreeVector( - OuterGeom->scintillator_y / 2., 0, 0), outer_steel_log, "HcalSteel", sandwich_log, 0, false, overlapcheck);
  G4RotationMatrix *hcal_rotm_sc = new G4RotationMatrix;
  hcal_rotm_sc->rotateZ(2 * M_PI*rad - OuterGeom->tilt_angle);
  outer_scinti = new G4PVPlacement(hcal_rotm_sc, G4ThreeVector(((OuterGeom->steel_plate_y_out + OuterGeom->steel_plate_y_in)/2.)/2.+0.0011*mm, 0, 0), outer_scinti_log, "HcalScinti", sandwich_log, 0, false, overlapcheck);
  return 0;
}

int
PHG4HcalTestBeamDetector::ConstructInnerSandwichVolume(G4LogicalVolume* sandwich_log)
{
  vector<G4LogicalVolume*> block_logic;
  G4VSolid*  inner_steel_trd = new G4Trap("Hcal_inner_steel_trd", inner_steel_z, inner_steel_x, inner_steel_y_out, inner_steel_y_in);
  G4LogicalVolume* inner_steel_log =  new G4LogicalVolume(inner_steel_trd, G4Material::GetMaterial("G4_Fe"), "Hcal_Inner_Steel", 0, 0, 0);
  G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetVisibility(false);
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Magenta());
  inner_steel_log->SetVisAttributes(cemcVisAtt);
  G4VSolid *inner_scintillator = new G4Box("Inner_Scinti_Box", inner_sc_dimension[1] / 2., inner_sc_dimension[0] / 2., inner_sc_dimension[2] / 2.);
  G4LogicalVolume* inner_scinti_log = new G4LogicalVolume(inner_scintillator, G4Material::GetMaterial("G4_POLYSTYRENE"), "Hcal_Inner_Sc", 0, 0, 0);
  G4VisAttributes* cemcVisAtt1 = new G4VisAttributes();
  cemcVisAtt1->SetVisibility(true);
  cemcVisAtt1->SetForceSolid(false);
  cemcVisAtt1->SetColour(G4Colour::Green());
  inner_scinti_log->SetVisAttributes(cemcVisAtt1);
  //  G4RotationMatrix *hcal_rotm = new G4RotationMatrix;
  //  hcal_rotm->rotateZ(270 * deg);
  // get to the bottom of the sandwich volome
  // -((y_out + y_in)/2)/2-sc/2
  // going 1/2 trapezoid up
  // -((y_out + y_in)/2)/2-sc/2 + ((y_out + y_in)/2)/2 = -sc/2

  inner_absorber = new G4PVPlacement(NULL, G4ThreeVector( - inner_sc_dimension[1] / 2., 0, 0), inner_steel_log, "HcalSteel", sandwich_log, 0, false, overlapcheck);
  G4RotationMatrix *hcal_rotm_sc = new G4RotationMatrix;
  hcal_rotm_sc->rotateZ((2 * M_PI - inner_tilt_angle)*rad);
  inner_scinti = new G4PVPlacement(hcal_rotm_sc, G4ThreeVector(((inner_steel_y_out + inner_steel_y_in)/2.)/2.+no_overlap, 0, 0), inner_scinti_log, "HcalScinti", sandwich_log, 0, false, overlapcheck);
  return 0;
}

void
PHG4HcalTestBeamDetector::CalculateGeometry()
{
  // tilt angle outer steel
  outer_tilt_angle =  1.41*M_PI/180.;
//atan((outer_steel_y_out-outer_steel_y_in)/(outer_plate_x+1*inch)); // add back 1 inch
// the inner size of the steel plate of the outer hcal in the drawings
// include the part used for bolting it to the inner hcal
// so we need to calculated it
  outer_steel_y_in = outer_steel_y_out - tan( outer_tilt_angle)*outer_plate_x;
  cout << "outer tilt angle: " << outer_tilt_angle*180./M_PI << ", outer_steel_y_in: " << outer_steel_y_in/inch << endl;
  inner_tilt_angle = 1.41*M_PI/180.;
  inner_steel_y_out = inner_steel_y_in + tan(inner_tilt_angle)*inner_steel_x;
    //atan((inner_steel_y_out-inner_steel_y_in)/inner_steel_x);
    cout << "inner tilt angle: " << inner_tilt_angle*180./M_PI << ", yout: " << inner_steel_y_out/inch 
	 << ", yin: " << inner_steel_y_in/inch << endl;
    //    inner_hcal_inner_radius = (inner_steel_y_out+inner_sc_dimension[1])/tan(inner_tilt_angle);
    //    outer_hcal_inner_radius = (outer_steel_y_in+outer_sc_dimension[1])/tan(outer_tilt_angle);
    inner_hcal_inner_radius = 188*cm;
    outer_hcal_inner_radius = 218*cm;
cout << "inner radius inner hcal: " << inner_hcal_inner_radius/cm << ", outer: " << outer_hcal_inner_radius/cm << endl;
    return;
}

HcalGeom::HcalGeom(const std::string &nam):
  name(nam)
{
  scinti_gap = 8.5*mm;
  if (name == "INNER")
    {
      InitInner();
    }
  else if (name == "OUTER")
    {
      InitOuter();
    }
  else
    {
      cout << "invalid name " << name << endl;
      exit(1);
    }
  CalculateGeometry();
  return;
}

void
HcalGeom::InitInner()
{
  // convert everything to G4 internal units
  num_sandwiches = 15;
  //  steel_plate_x = 12.43*inch;
  steel_plate_x = 13.89*inch;
  steel_plate_y_in = 0.787*inch;
  steel_plate_y_out = 1.115*inch;
  steel_plate_z = 17.1*inch;
  scintillator_x = 12.43*inch-10*mm;
  scintillator_y =  7*mm;
  scintillator_z = steel_plate_z;
  tilt_angle =  1.41*deg;
  CalculateGeometryInner();
}

void
HcalGeom::InitOuter()
{
  num_sandwiches = 15;
  steel_plate_x = 27.1 * inch;
  steel_plate_y_in = 1.099 * inch;
  steel_plate_y_out = 1.776 * inch;
  steel_plate_z = 29.53 * inch;
  scintillator_x =  649.8*mm;
  scintillator_y =  7 * mm;
  scintillator_z = steel_plate_z ;
  CalculateGeometryOuter();
}

void
HcalGeom::CalculateGeometryInner()
{
  // inner hcal has angle (1.41 deg but not steel_plate_y_out)
  steel_plate_y_out = (tan(tilt_angle*rad)*steel_plate_x)+ steel_plate_y_in;
 cout << name << ": tilt angle: " << tilt_angle/deg << endl;
 cout << name << ": steel_plate_y_out: " << steel_plate_y_out/inch << endl;
 cout << name << ": inner radius: " << inner_radius/cm << endl;
    return;
}

void
HcalGeom::CalculateGeometryOuter()
{
  // tilt angle outer steel
  //  tilt_angle =  1.41*M_PI/180.;
tilt_angle = atan((steel_plate_y_out-steel_plate_y_in)/steel_plate_x)*rad;
 cout << name << ": tilt angle: " << tilt_angle/deg << endl;
    return;
}

void
HcalGeom::CalculateGeometry()
{
  scintillator_y_vertical_projection  = fabs(scintillator_y/cos(M_PI-tilt_angle/rad));
  double steelscinti_y = scintillator_y_vertical_projection + steel_plate_y_in;

  inner_radius = steelscinti_y/tan(tilt_angle/rad);
 cout << name << ": steelscinti_y: " << steelscinti_y/inch << endl;
 cout << name << ": inner radius: " << inner_radius/cm << endl;
  // tilt angle outer steel
  //  tilt_angle =  1.41*M_PI/180.;
 cout << name << ": tilt angle: " << tilt_angle/deg << endl;
 cout << name << ": inner radius: " << inner_radius/cm << endl;
    return;
}

int
PHG4HcalTestBeamDetector::DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol,G4RotationMatrix *rotm )
{
  static int i = 0;
  G4LogicalVolume* checksolid = new G4LogicalVolume(volume,G4Material::GetMaterial("G4_POLYSTYRENE"),"HcalOuterScinti", 0, 0, 0);
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
    default:
      visattchk->SetColour(G4Colour::Green());
      break;
    }

  checksolid->SetVisAttributes(visattchk);
  new G4PVPlacement(rotm,G4ThreeVector(0,0,0),checksolid,"CHKVOL",logvol, 0, false, overlapcheck);
  return 0;
}
