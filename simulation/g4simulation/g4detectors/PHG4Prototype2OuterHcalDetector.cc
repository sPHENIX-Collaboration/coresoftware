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
#include <CGAL/Boolean_set_operations_2.h>

#include <boost/math/special_functions/sign.hpp>

#include <cmath>
#include <sstream>

typedef CGAL::Exact_circular_kernel_2             Circular_k;
typedef CGAL::Point_2<Circular_k>                 Point_2;
typedef CGAL::Circle_2<Circular_k>                Circle_2;
typedef CGAL::Circular_arc_point_2<Circular_k>          Circular_arc_point_2;
typedef CGAL::Line_2<Circular_k>                Line_2;
typedef CGAL::Segment_2<Circular_k>                Segment_2;

using namespace std;

// there is still a minute problem for very low tilt angles where the scintillator
// face touches the boundary instead of the corner, subtracting 1 permille from the total
// scintilator length takes care of this
static double subtract_from_scinti_x = 0.1*mm;

PHG4Prototype2OuterHcalDetector::PHG4Prototype2OuterHcalDetector( PHCompositeNode *Node, PHG4Parameters *parameters, const std::string &dnam  ):
  PHG4Detector(Node, dnam),
  params(parameters),
  outerhcalsteelplate(NULL),
  inner_radius(1830*mm),
  outer_radius(2685*mm),
  steel_x(823.*mm),
  steel_yhi(42.5*mm),
  steel_ylo(26.2*mm),
  steel_z(1600.*mm),
  bottom_xmiddle_steel_tile((2601.2*mm-1777.6*mm)/2.+1777.6*mm),
  //  bottom_ymiddle_steel_tile(-459.8*mm+(steel_yhi+steel_ylo)/4.),
  bottom_ymiddle_steel_tile(-459.8*mm),
  size_z(1600*mm),
  scinti_tile_x(NAN),
  scinti_tile_x_lower(NAN),
  scinti_tile_x_upper(NAN),
  scinti_tile_z(size_z),
  scinti_tile_thickness(params->get_double_param("scinti_tile_thickness")*cm),
  scinti_gap(8.482*mm),
  tilt_angle(12*deg),
  envelope_inner_radius(inner_radius),
  envelope_outer_radius(outer_radius),
  envelope_z(size_z),
  volume_envelope(NAN),
  volume_steel(NAN),
  volume_scintillator(NAN),
  n_scinti_plates(20),
  n_scinti_tiles(params->get_int_param("n_scinti_tiles")),
  active(params->get_int_param("active")),
  absorberactive(params->get_int_param("absorberactive")),
  layer(0),
  scintilogicnameprefix("HcalInnerScinti")
{

  // allocate memory for scintillator plates
  scinti_tiles_vec.assign(2 * n_scinti_tiles, static_cast<G4VSolid *>(NULL));
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
      if (steel_absorber_vec.find(volume) != steel_absorber_vec.end())
	{
	  return -1;
	}
    }
  if (active)
    {
      if (volume->GetName().find(scintilogicnameprefix) != string::npos)
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
      G4VSolid* steel_plate = new G4Trap("OuterHcalSteelPlate",steel_z,steel_x,steel_yhi,steel_ylo);
      volume_steel = steel_plate->GetCubicVolume()*n_scinti_plates;
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
PHG4Prototype2OuterHcalDetector::ConstructSteelScintiVolume(G4LogicalVolume* hcalenvelope)
{
  G4LogicalVolume* steel_plate = ConstructSteelPlate(hcalenvelope);
  G4VisAttributes* visattchk = new G4VisAttributes();

  G4LogicalVolume* scintibox_logical = ConstructScintillatorBox(hcalenvelope);
  double yhi = steel_yhi+scinti_gap;
  double ylo = steel_ylo+scinti_gap;

  G4VSolid* steelscintimother = new  G4Trap("steelscintimother",steel_z,steel_x,yhi,ylo);
  G4LogicalVolume* steelscintilogical = new G4LogicalVolume(steelscintimother,G4Material::GetMaterial("G4_AIR"),G4String("steelscintimotherlog"), 0, 0, 0);
  G4RotationMatrix *rot = new G4RotationMatrix();
  double ymiddle = scinti_gap + (steel_yhi + steel_ylo)/2.;
  steel_absorber_vec.insert(new G4PVPlacement(rot, G4ThreeVector(-(yhi+ylo)/4.+scinti_gap+(steel_yhi + steel_ylo)/4.,0,0),steel_plate,"steel",steelscintilogical,false,0,overlapcheck));
  double ydown = ymiddle; 
  rot = new G4RotationMatrix();
  new G4PVPlacement(rot, G4ThreeVector(-(yhi+ylo)/4.+scinti_gap/2.,0,0),scintibox_logical,"scintibox",steelscintilogical,false,0,overlapcheck);
  visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  visattchk->SetColour(G4Colour::Green());
  steelscintilogical->SetVisAttributes(visattchk);
  return steelscintilogical;


}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintillatorBox(G4LogicalVolume* hcalenvelope)
{ 
  G4VSolid* scintiboxsolid = new G4Box("ScintiBox_Solid",scinti_gap/2.,steel_x/2.,scinti_tile_z/2.);
  //  DisplayVolume(scintiboxsolid,hcalenvelope);
  G4LogicalVolume* scintiboxlogical = new G4LogicalVolume(scintiboxsolid,G4Material::GetMaterial("G4_AIR"),G4String("Hcal_scintibox"), 0, 0, 0);
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  visattchk->SetColour(G4Colour::Yellow());
  scintiboxlogical->SetVisAttributes(visattchk);
  ConstructScintiTile_1(hcalenvelope);
  return scintiboxlogical;
}

G4VSolid*
PHG4Prototype2OuterHcalDetector::ConstructScintiTile_1(G4LogicalVolume* hcalenvelope)
{
  G4double outer2UpDz = 0.5*649.8*mm;
  G4double outer2UpTheta = 8.8*M_PI/180.;
  G4double outer2UpPhi = 0.0*M_PI/180.;
  G4double outer2UpDy1 = 0.35*cm;
  G4double outer2UpDy2 = 0.35*cm;
  G4double outer2UpDx1 = 0.5*179.3*mm, outer2UpDx2 = 0.5*179.3*mm;
  G4double outer2UpDx3 = 0.5*113.2*mm, outer2UpDx4 = 0.5*113.2*mm;
  G4double outer2UpAlp1 = 0.*M_PI/180., outer2UpAlp2 = 0.*M_PI/180;

  G4VSolid *outer2USheetSolid = new G4Trap("outer2USheet",
				  outer2UpDz,
				  outer2UpTheta,
				  outer2UpPhi, 
				  outer2UpDy1,
				  outer2UpDx1,
				  outer2UpDx2,
				  outer2UpAlp1,
				  outer2UpDy2,
				  outer2UpDx3,
				  outer2UpDx4,
				  outer2UpAlp2);
  DisplayVolume(outer2USheetSolid,hcalenvelope);
  return outer2USheetSolid;

  double xbase = 1770.9;
  double ybase = -460.3;
  G4TwoVector v4(xbase*mm,ybase*mm);
  G4TwoVector v3((xbase+828.9)*mm,ybase*mm);
  G4TwoVector v2((xbase+828.9)*mm,(ybase-240.54)*mm);
  G4TwoVector v1(xbase*mm,(ybase-166.2)*mm);
  cout << "v1: " << v1 << endl;
  cout << "v2: " << v2 << endl;
  cout << "v3: " << v3 << endl;
  cout << "v4: " << v4 << endl;
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(v1);
  vertexes.push_back(v2);
  vertexes.push_back(v3);
  vertexes.push_back(v4);
  G4TwoVector zero(0, 0);
  G4VSolid* steel_plate =  new G4ExtrudedSolid("ScintiTile_1",
					       vertexes,
					       7*mm  / 2.0,
					       zero, 1.0,
					       zero, 1.0);

  //  G4RotationMatrix *rotm = new G4RotationMatrix();
  // rotm->rotateX(90*deg);
  // DisplayVolume(steel_plate, hcalenvelope, rotm);
  //volume_steel = steel_plate->GetCubicVolume()*n_scinti_plates;
  return steel_plate;
}

// Construct the envelope and the call the
// actual inner hcal construction
void
PHG4Prototype2OuterHcalDetector::Construct( G4LogicalVolume* logicWorld )
{
  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  //  G4VSolid* hcal_envelope_cylinder = new G4Tubs("OuterHcal_envelope_solid",  envelope_inner_radius-0.5*cm, envelope_outer_radius, envelope_z / 2., tan(-459.8/1770.9), fabs(tan(-459.8/1770.9)) +  12.4/180.* M_PI);
  G4VSolid* hcal_envelope_cylinder = new G4Tubs("OuterHcal_envelope_solid",  0*cm, envelope_outer_radius, envelope_z / 2., 0, 2.* M_PI);
  volume_envelope = hcal_envelope_cylinder->GetCubicVolume();
  G4LogicalVolume* hcal_envelope_log =  new G4LogicalVolume(hcal_envelope_cylinder, Air, G4String("Hcal_envelope"), 0, 0, 0);
  G4VisAttributes* hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::White());
  hcal_envelope_log->SetVisAttributes(hcalVisAtt);
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(params->get_double_param("rot_x")*deg);
  hcal_rotm.rotateY(params->get_double_param("rot_y")*deg);
  hcal_rotm.rotateZ(params->get_double_param("rot_z")*deg);
  //  new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(0,0,0)), hcal_envelope_log, "OuterHcalEnvelope", logicWorld, 0, false, overlapcheck);
  new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(params->get_double_param("place_x")*cm, params->get_double_param("place_y")*cm, params->get_double_param("place_z")*cm)), hcal_envelope_log, "OuterHcalEnvelope", logicWorld, 0, false, overlapcheck);
  //ConstructSteelScintiVolume(hcal_envelope_log);
  //ConstructSteelPlate(hcal_envelope_log);
  //ConstructScintiTile_1(hcal_envelope_log);
  ConstructOuterHcal(hcal_envelope_log);
  //  AddGeometryNode();
  return;
}

int
PHG4Prototype2OuterHcalDetector::ConstructOuterHcal(G4LogicalVolume* hcalenvelope)
{
  G4LogicalVolume* steel_plate = ConstructSteelPlate(hcalenvelope); // bottom steel plate
  G4VisAttributes* visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  visattchk->SetColour(G4Colour::Blue());
  steel_plate->SetVisAttributes(visattchk);
  G4LogicalVolume  *steelscinti = ConstructSteelScintiVolume(hcalenvelope);
  visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  visattchk->SetColour(G4Colour::Grey());
  steelscinti->SetVisAttributes(visattchk);
  double phi = 0.;
  double phislat = 0.;
  double deltaphi = 2 * M_PI / n_scinti_plates;
  deltaphi = 2 * M_PI / 320.;
  ostringstream name;
  //  double middlerad = outer_radius - (outer_radius - inner_radius) / 2.;
  double shiftslat = fabs(scinti_tile_x_lower - scinti_tile_x_upper)/2.;
  double bottomslat_y = -460.3*mm;
  name.str("");
  name << "Steel_plate_1";
  G4RotationMatrix *Rot = new G4RotationMatrix();
  Rot->rotateZ(-90*deg);
  steel_absorber_vec.insert(new G4PVPlacement(Rot, G4ThreeVector(bottom_xmiddle_steel_tile, bottom_ymiddle_steel_tile, 0), steel_plate, name.str().c_str(), hcalenvelope, false, 0, overlapcheck));
  return 0;
  //  double middlerad = sqrt(bottom_xmiddle_steel_tile*bottom_xmiddle_steel_tile + (bottom_ymiddle_steel_tile+scinti_gap/2.) * (bottom_ymiddle_steel_tile+scinti_gap/2.));
  double middlerad = sqrt(bottom_xmiddle_steel_tile*bottom_xmiddle_steel_tile + (bottom_ymiddle_steel_tile) * (bottom_ymiddle_steel_tile));
  //  double philow = atan((bottom_ymiddle_steel_tile-scinti_gap/2.)/bottom_xmiddle_steel_tile);
  double philow = atan((bottom_ymiddle_steel_tile)/bottom_xmiddle_steel_tile);
  double steelrat = ((steel_yhi + steel_ylo)/2.)/((steel_yhi + steel_ylo)/2. + scinti_gap);
  cout << "steel ratio: " << steelrat << " deltaphi: " << deltaphi*180/M_PI << ", steel phi: " << deltaphi*steelrat*180./M_PI << endl;
  //  phi -= deltaphi*steelrat;
  cout << "philow: " << philow*180./M_PI << endl;
  //return 0;
  //phi -= deltaphi;
  //  phi+= deltaphi;//*steelrat;
  phi+= 1.032*M_PI/180.;//*steelrat;
  for (int i = 0; i < n_scinti_plates; i++)
    //for (int i = 1; i < 2; i++)
    {
      double ypos = sin(phi+philow) * middlerad;
      double xpos = cos(phi+philow) * middlerad;
      cout << "bottom_ymiddle_steel_tile: " << bottom_ymiddle_steel_tile
	   << ", ypos: " << ypos << endl;
      cout << "bottom_xmiddle_steel_tile: " << bottom_xmiddle_steel_tile
	   << ", xpos: " << xpos << endl;
      // the center of the scintillator is not the center of the inner hcal
      // but depends on the tilt angle. Therefore we need to shift
      // the center from the mid point
      // ypos += sin((-tilt_angle)/rad - phi)*shiftslat;
      // xpos -= cos((-tilt_angle)/rad - phi)*shiftslat;
      Rot = new G4RotationMatrix();
      Rot->rotateZ(-90*deg - phi);
            G4ThreeVector g4vec(xpos, ypos, 0);
	    //G4ThreeVector g4vec(bottom_xmiddle_steel_tile, bottom_ymiddle_steel_tile, 0);
      //      scinti_mother_logical->MakeImprint(hcalenvelope, g4vec, Rot, i, overlapcheck);
      name.str("");
      name << "OuterHcalSteel_" << i;
      steel_absorber_vec.insert(new G4PVPlacement(Rot, g4vec, steelscinti, name.str().c_str(), hcalenvelope, false, i, overlapcheck));
      bottomslat_y += 30*mm;
      phislat += deltaphi;
      phi += deltaphi;
    }
  return 0;
}


// split the big scintillator into tiles covering eta ranges
// since they are tilted it is not a straightforward theta cut
// it starts at a given eta at the inner radius but the outer radius needs adjusting
void
PHG4Prototype2OuterHcalDetector::ConstructHcalSingleScintillators(G4LogicalVolume* hcalenvelope)
{
  return;
}

G4AssemblyVolume *
PHG4Prototype2OuterHcalDetector::ConstructHcalScintillatorAssembly(G4LogicalVolume* hcalenvelope)
{
  ConstructHcalSingleScintillators(hcalenvelope);
  G4AssemblyVolume *assmeblyvol = new G4AssemblyVolume();
  ostringstream name;
  G4ThreeVector g4vec;

  double steplimits = params->get_double_param("steplimits")*cm;
  for (unsigned int i = 0; i < scinti_tiles_vec.size(); i++)
    {
      name.str("");
      name << scintilogicnameprefix << i;
      G4UserLimits *g4userlimits = NULL;
      if (isfinite(steplimits))
	{
	  g4userlimits = new G4UserLimits(steplimits);
	}
      G4LogicalVolume *scinti_tile_logic = new G4LogicalVolume(scinti_tiles_vec[i], G4Material::GetMaterial("G4_POLYSTYRENE"), name.str().c_str(), NULL, NULL, g4userlimits);
      G4VisAttributes *visattchk = new G4VisAttributes();
      visattchk->SetVisibility(true);
      visattchk->SetForceSolid(true);
      visattchk->SetColour(G4Colour::Green());
      scinti_tile_logic->SetVisAttributes(visattchk);
      assmeblyvol->AddPlacedVolume(scinti_tile_logic, g4vec, NULL);
    }
  return assmeblyvol;
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
  cout << "Inner Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
    {
      cout << "Volume Envelope: " << volume_envelope/cm/cm/cm << " cm^3" << endl;
      cout << "Volume Steel: " << volume_steel/cm/cm/cm << " cm^3" << endl;
      cout << "Volume Scintillator: " << volume_scintillator/cm/cm/cm << " cm^3" << endl;
      cout << "Volume Air: " << (volume_envelope - volume_steel - volume_scintillator)/cm/cm/cm << " cm^3" << endl;
    }
  return;
}
