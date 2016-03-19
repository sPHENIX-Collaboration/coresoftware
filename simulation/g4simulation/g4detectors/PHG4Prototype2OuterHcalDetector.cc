#include "PHG4Prototype2OuterHcalDetector.h"
#include "PHG4Parameters.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeomv3.h"

#include <g4main/PHG4Utils.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Object.h>
#include <CGAL/Circular_kernel_intersections.h>

#include <boost/math/special_functions/sign.hpp>

#include <cmath>
#include <sstream>

typedef CGAL::Exact_circular_kernel_2             Circular_k;
typedef CGAL::Point_2<Circular_k>                 Point_2;
typedef CGAL::Line_2<Circular_k>                Line_2;

using namespace std;

PHG4Prototype2OuterHcalDetector::PHG4Prototype2OuterHcalDetector( PHCompositeNode *Node, PHG4Parameters *parameters, const std::string &dnam  ):
  PHG4Detector(Node, dnam),
  params(parameters),
  outerhcalsteelplate(NULL),
  outerhcalassembly(NULL),
  steel_plate_corner_upper_left(1777.6*mm,-433.5*mm),
  steel_plate_corner_upper_right(2600.4*mm,-417.4*mm), 
  steel_plate_corner_lower_right(2601.2*mm,-459.8*mm),
  steel_plate_corner_lower_left(1770.9*mm,-459.8*mm),
  scinti_u1_front_size(166.2*mm),
  scinti_u1_corner_upper_left(0*mm,0*mm),
  scinti_u1_corner_upper_right(828.9*mm,0*mm),
  scinti_u1_corner_lower_right(828.9*mm,-240.54*mm),
  scinti_u1_corner_lower_left(0*mm,-scinti_u1_front_size),
  scinti_u2_corner_upper_left(0*mm,0*mm),
  scinti_u2_corner_upper_right(828.9*mm,-74.3*mm),
  scinti_u2_corner_lower_right(828.9*mm,-320.44*mm),
  scinti_u2_corner_lower_left(0*mm,-171.0*mm),
  inner_radius(1830*mm),
  outer_radius(2685*mm),
  scinti_x(828.9),
  steel_x(823.*mm),
  steel_z(1600.*mm),
  size_z(steel_z),
  scinti_tile_z(steel_z),
  scinti_tile_thickness(7*mm),
  scinti_box_shift(1.09*mm), // that was found experimetnally by removing overlaps
  gap_between_tiles(1*mm),
  scinti_gap(8.5*mm),
  tilt_angle(12*deg),
  deltaphi(2*M_PI/320.),
  volume_steel(NAN),
  volume_scintillator(NAN),
  n_scinti_plates(20),
  n_steel_plates(n_scinti_plates+1),
  active(params->get_int_param("active")),
  absorberactive(params->get_int_param("absorberactive")),
  layer(0),
  scintilogicnameprefix("OuterHcalScintiMother")
{
}

//_______________________________________________________________
//_______________________________________________________________
int
PHG4Prototype2OuterHcalDetector::IsInPrototype2OuterHcal(G4VPhysicalVolume * volume) const
{
  // G4AssemblyVolumes naming convention:
  //     av_WWW_impr_XXX_YYY_ZZZ

  // where:

  //     WWW - assembly volume instance number
  //     XXX - assembly volume imprint number
  //     YYY - the name of the placed logical volume
  //     ZZZ - the logical volume index inside the assembly volume
  // e.g. av_1_impr_82_HcalInnerScinti_11_pv_11
  // 82 the number of the scintillator mother volume
  // HcalInnerScinti_11: name of scintillator slat
  // 11: number of scintillator slat logical volume
  if (absorberactive)
    {
      if (volume->GetName().find("OuterHcalSteelPlate") != string::npos)
	{
	  return -1;
	}
    }
  if (active)
    {
      if (volume->GetName().find("OuterScinti") != string::npos)
	{
	  return 1;
	}
    }
  return 0;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructSteelPlate(G4LogicalVolume* hcalenvelope)
{
  if (!outerhcalsteelplate)
    {
      G4VSolid* steel_plate;
      std::vector<G4TwoVector> vertexes;
      vertexes.push_back(steel_plate_corner_upper_left);
      vertexes.push_back(steel_plate_corner_upper_right);
      vertexes.push_back(steel_plate_corner_lower_right);
      vertexes.push_back(steel_plate_corner_lower_left);
      G4TwoVector zero(0, 0);
      steel_plate =  new G4ExtrudedSolid("OuterHcalSteelPlateSolid",
					 vertexes,
					 size_z  / 2.0,
					 zero, 1.0,
					 zero, 1.0);

      volume_steel = steel_plate->GetCubicVolume()*n_steel_plates;
      outerhcalsteelplate = new G4LogicalVolume(steel_plate,G4Material::GetMaterial("SS310"),"OuterHcalSteelPlate", 0, 0, 0);
      G4VisAttributes* visattchk = new G4VisAttributes();
      visattchk->SetVisibility(true);
      visattchk->SetForceSolid(false);
      visattchk->SetColour(G4Colour::Blue());
      outerhcalsteelplate->SetVisAttributes(visattchk);
    }
  return outerhcalsteelplate;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintillatorBox(G4LogicalVolume* hcalenvelope)
{ 
  G4VSolid* scintiboxsolid = new G4Box("OuterHcalScintiMother",scinti_x/2.,scinti_gap/2.,scinti_tile_z/2.);
  //  DisplayVolume(scintiboxsolid,hcalenvelope);
  G4LogicalVolume* scintiboxlogical = new G4LogicalVolume(scintiboxsolid,G4Material::GetMaterial("G4_AIR"),G4String("OuterHcalScintiMother"), 0, 0, 0);
  G4VisAttributes* hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Red());
  G4LogicalVolume *scintiu1_logic = ConstructScintiTileU1(hcalenvelope);
  scintiu1_logic->SetVisAttributes(hcalVisAtt);

  hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Cyan());
  G4LogicalVolume *scintiu2_logic = ConstructScintiTileU2(hcalenvelope);
  scintiu2_logic->SetVisAttributes(hcalVisAtt);
  G4RotationMatrix *Rot;  
  Rot = new G4RotationMatrix();  
  Rot->rotateX(-90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,-scinti_u1_front_size-gap_between_tiles/2.-gap_between_tiles),scintiu2_logic,"OuterScinti_0", scintiboxlogical, false, 0, overlapcheck);

  Rot = new G4RotationMatrix();  
  Rot->rotateX(-90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,-gap_between_tiles/2.),scintiu1_logic,"OuterScinti_1", scintiboxlogical, false, 0, overlapcheck);

  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,gap_between_tiles/2.),scintiu1_logic,"OuterScinti_2", scintiboxlogical, false, 0, overlapcheck);

  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,scinti_u1_front_size+gap_between_tiles/2.+gap_between_tiles),scintiu2_logic,"OuterScinti_3", scintiboxlogical, false, 0, overlapcheck);


  return scintiboxlogical;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTileU1(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_u1_corner_upper_left);
  vertexes.push_back(scinti_u1_corner_upper_right);
  vertexes.push_back(scinti_u1_corner_lower_right);
  vertexes.push_back(scinti_u1_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintiu1 =  new G4ExtrudedSolid("OuterHcalScintiU1",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintiu1_logic = new G4LogicalVolume(scintiu1,G4Material::GetMaterial("G4_POLYSTYRENE"),"OuterHcalScintiU1", NULL, NULL, NULL);
  //   DisplayVolume(scintiu1,hcalenvelope);
  return scintiu1_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTileU2(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_u2_corner_upper_left);
  vertexes.push_back(scinti_u2_corner_upper_right);
  vertexes.push_back(scinti_u2_corner_lower_right);
  vertexes.push_back(scinti_u2_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintiu2 =  new G4ExtrudedSolid("OuterHcalScintiU2",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintiu2_logic = new G4LogicalVolume(scintiu2,G4Material::GetMaterial("G4_POLYSTYRENE"),"OuterHcalScintiU2", NULL, NULL, NULL);
  //   DisplayVolume(scintiu2,hcalenvelope);
  return scintiu2_logic;
}

// Construct the envelope and the call the
// actual inner hcal construction
void
PHG4Prototype2OuterHcalDetector::Construct( G4LogicalVolume* logicWorld )
{
  G4ThreeVector g4vec(0,0,0);
  G4RotationMatrix *Rot = new G4RotationMatrix();
  Rot->rotateX(params->get_double_param("rot_x")*deg);
  Rot->rotateY(params->get_double_param("rot_y")*deg);
  Rot->rotateZ(params->get_double_param("rot_z")*deg);
  outerhcalassembly = new G4AssemblyVolume();
  //ConstructSteelPlate(hcal_envelope_log);
  // return;
  ConstructOuterHcal(logicWorld);
  outerhcalassembly->MakeImprint(logicWorld,g4vec,Rot,0,overlapcheck);
  //  AddGeometryNode();
  return;
}

int
PHG4Prototype2OuterHcalDetector::ConstructOuterHcal(G4LogicalVolume* hcalenvelope)
{
  G4LogicalVolume* steel_plate = ConstructSteelPlate(hcalenvelope); // bottom steel plate
  G4LogicalVolume* scintibox = ConstructScintillatorBox(hcalenvelope);
  double phi = 0.;
  double phislat = 0.;
  ostringstream name;
  // the coordinate of the center of the bottom of the bottom steel plate
  // to get the radius of the circle which is the center of the scintillator box
  double bottom_xmiddle_steel_tile = (steel_plate_corner_lower_right.x()-steel_plate_corner_lower_left.x())/2.+steel_plate_corner_lower_left.x();
  double bottom_ymiddle_steel_tile = steel_plate_corner_lower_right.y();
  double middlerad = sqrt(bottom_xmiddle_steel_tile*bottom_xmiddle_steel_tile + bottom_ymiddle_steel_tile * bottom_ymiddle_steel_tile);
  double philow = atan((bottom_ymiddle_steel_tile-scinti_gap/2.)/bottom_xmiddle_steel_tile);
  double scintiangle = GetScintiAngle();
  for (int i = 0; i < n_steel_plates; i++)
    //      for (int i = 0; i < 2; i++)
    {
      name.str("");
      name << "OuterHcalSteel_" << i;
      G4RotationMatrix *Rot = new G4RotationMatrix();
      Rot->rotateZ(phi*rad);
      G4ThreeVector g4vec(0,0,0);
      outerhcalassembly->AddPlacedVolume(steel_plate,g4vec,Rot);
      if (i > 0)
	{
	  double ypos = sin(phi+philow) * middlerad;
	  double xpos = cos(phi+philow) * middlerad;
	  // the center of the scintillator is not the center of the inner hcal
	  // but depends on the tilt angle. Therefore we need to shift
	  // the center from the mid point
	  ypos += sin((-tilt_angle)/rad - phi)*scinti_box_shift;
	  xpos -= cos((-tilt_angle)/rad - phi)*scinti_box_shift;
	  name.str("");
	  name << "OuterHcalScintiBox_" << i;
	  Rot = new G4RotationMatrix();
	  Rot->rotateZ(scintiangle+phislat);
	  G4ThreeVector g4vec(xpos, ypos, 0);

	  outerhcalassembly->AddPlacedVolume(scintibox,g4vec,Rot);
	  phislat += deltaphi;
	}
      phi += deltaphi;
    }
  return 0;
}

// calculate the angle of the bottom scintillator. It is the angle of the top edge
// of the steel plate
double
PHG4Prototype2OuterHcalDetector::GetScintiAngle()
{
  Point_2 upleft(steel_plate_corner_upper_left.x(),steel_plate_corner_upper_left.y());
  // we just need a horizontal line from the upper left corner for the intersection
  // with the rear end of the steel plate
  Point_2 lefttmp(steel_plate_corner_upper_left.x()+200*mm,steel_plate_corner_upper_left.y());
  Point_2 upright(steel_plate_corner_upper_right.x(),steel_plate_corner_upper_right.y());
  Point_2 downright(steel_plate_corner_lower_right.x(),steel_plate_corner_lower_right.y());
  Line_2 rightside(upright,downright);
  Line_2 horiz(upleft,lefttmp);
  CGAL::Object result = CGAL::intersection(rightside, horiz);
  // this point is the intersection of the horizontal line going through the
  // upper left steel plate cornet with the rear edge of the steel plate
  Point_2 intersect;
  if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
    {
      intersect = *ipoint;
    }
  // length of rear edge of steel plate from upper right corner to intersection point
  double lenshort = sqrt((steel_plate_corner_upper_right.x()-CGAL::to_double(intersect.x()))
                        *(steel_plate_corner_upper_right.x()-CGAL::to_double(intersect.x()))
		       + (steel_plate_corner_upper_right.y()-CGAL::to_double(intersect.y()))
                        *(steel_plate_corner_upper_right.y()-CGAL::to_double(intersect.y())));
  double lenup = sqrt((steel_plate_corner_upper_right.x() - steel_plate_corner_upper_left.x())
                     *(steel_plate_corner_upper_right.x() - steel_plate_corner_upper_left.x())
		    + (steel_plate_corner_upper_right.y() - steel_plate_corner_upper_left.y())
		      *(steel_plate_corner_upper_right.y() - steel_plate_corner_upper_left.y()));
  double angle = asin(lenshort/lenup);
  return angle;
}

int
PHG4Prototype2OuterHcalDetector::DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix *rotm )
{
  static int i = 0;
  G4LogicalVolume* checksolid = new G4LogicalVolume(volume, G4Material::GetMaterial("G4_POLYSTYRENE"), "DISPLAYLOGICAL", 0, 0, 0);
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
      i = 0;
      break;
    }

  checksolid->SetVisAttributes(visattchk);
  new G4PVPlacement(rotm, G4ThreeVector(0, 0, 0), checksolid, "DISPLAYVOL", logvol, 0, false, overlapcheck);
  //  new G4PVPlacement(rotm, G4ThreeVector(0, -460.3, 0), checksolid, "DISPLAYVOL", logvol, 0, false, overlapcheck);
  return 0;
}

void
PHG4Prototype2OuterHcalDetector::AddGeometryNode()
{
  if (params->get_int_param("active"))
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
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeomv3(inner_radius / cm, (params->get_double_param("place_z")*cm - size_z / 2.) / cm, (params->get_double_param("place_z")*cm + size_z / 2.) / cm, (outer_radius - inner_radius) / cm, n_scinti_plates,  tilt_angle / rad, 0);
      geo->AddLayerGeom(layer, mygeom);
      if (verbosity > 0) geo->identify();
    }
}

void
PHG4Prototype2OuterHcalDetector::Print(const string &what) const
{
  cout << "Outer Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
    {
      cout << "Volume Steel: " << volume_steel/cm/cm/cm << " cm^3" << endl;
      cout << "Volume Scintillator: " << volume_scintillator/cm/cm/cm << " cm^3" << endl;
    }
  return;
}
