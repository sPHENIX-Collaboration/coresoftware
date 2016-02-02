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
  inner_radius(1830*mm),
  outer_radius(2685*mm),
  size_z(1600*mm),
  scinti_tile_x(NAN),
  scinti_tile_x_lower(NAN),
  scinti_tile_x_upper(NAN),
  scinti_tile_z(size_z),
  scinti_tile_thickness(params->get_double_param("scinti_tile_thickness")*cm),
  scinti_gap(params->get_double_param("scinti_gap")*cm),
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

G4VSolid*
PHG4Prototype2OuterHcalDetector::ConstructScintillatorBox(G4LogicalVolume* hcalenvelope)
{
  double mid_radius = inner_radius + (outer_radius - inner_radius) / 2.;
  Point_2 p_in_1(mid_radius, 0); // center of scintillator

  // length of upper edge (middle till outer circle intersect
  // x/y coordinate of end of center vertical
  double xcoord  = scinti_tile_thickness / 2. *sin(fabs(tilt_angle)/rad) +  mid_radius;
  double ycoord  =   scinti_tile_thickness / 2. * cos(fabs(tilt_angle)/rad) + 0;
  Point_2 p_upperedge(xcoord, ycoord);
  Line_2 s2(p_in_1, p_upperedge); // center vertical

  Line_2 perp =  s2.perpendicular(p_upperedge);
  Point_2 sc1(outer_radius, 0), sc2(0, outer_radius), sc3(-outer_radius, 0);
  Circle_2 outer_circle(sc1, sc2, sc3);
  vector< CGAL::Object > res;
  CGAL::intersection(outer_circle, perp, std::back_inserter(res));
  Point_2 upperright;
  vector< CGAL::Object >::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); ++iter)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) >  CGAL::to_double(p_upperedge.x()))
	    {
	      double deltax = CGAL::to_double(point->first.x()) - CGAL::to_double(p_upperedge.x());
	      double deltay = CGAL::to_double(point->first.y()) - CGAL::to_double(p_upperedge.y());
	      // the scintillator is twice as long
	      scinti_tile_x_upper = sqrt(deltax * deltax + deltay * deltay); //
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      upperright = pntmp;
	    }
	}
      else
	{
	  cout << "CGAL::Object type not pair..." << endl;
	}
    }
  // length of lower edge (middle till inner circle intersect
  xcoord = mid_radius - scinti_tile_thickness / 2. * sin(fabs(tilt_angle)/rad);
  ycoord  = 0 -  scinti_tile_thickness / 2. * cos(fabs(tilt_angle)/rad);
  Point_2 p_loweredge(xcoord, ycoord);
  Line_2 s3(p_in_1, p_loweredge);
  Line_2 l_lower = s3.perpendicular(p_loweredge);
  Point_2 ic1(inner_radius, 0), ic2(0, inner_radius), ic3(-inner_radius, 0);
  Circle_2 inner_circle(ic1,ic2,ic3);
  res.clear();
  CGAL::intersection(inner_circle,l_lower, std::back_inserter(res));
  Point_2 lowerleft;
  // we have 2 intersections - we want the one furthest to the right (largest x). The correct one is
  // certainly > 0 but the second one depends on the tilt angle and might also be > 0
  double minx = 0;
  for (iter = res.begin(); iter != res.end(); ++iter)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) >  minx)
	    {
	      minx = CGAL::to_double(point->first.x());
	      double deltax = CGAL::to_double(point->first.x()) - CGAL::to_double(p_loweredge.x());
	      double deltay = CGAL::to_double(point->first.y()) - CGAL::to_double(p_loweredge.y());
              scinti_tile_x_lower = sqrt(deltax * deltax + deltay * deltay);
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      lowerleft = pntmp;
	    }
	}
    }
  scinti_tile_x  = scinti_tile_x_upper + scinti_tile_x_lower;
  scinti_tile_x  -= subtract_from_scinti_x;
  G4VSolid* scintibox =  new G4Box("ScintiTile", scinti_tile_x / 2., scinti_tile_thickness / 2., scinti_tile_z / 2.);
  volume_scintillator = scintibox->GetCubicVolume() *n_scinti_plates;

  return scintibox;
}

G4VSolid*
PHG4Prototype2OuterHcalDetector::ConstructSteelPlate(G4LogicalVolume* hcalenvelope)
{
  // calculate steel plate on top of the scinti box. Lower edge is the upper edge of
  G4TwoVector v4(1770.9*mm,-459.8*mm);
  G4TwoVector v3(2601.2*mm,-459.8*mm);
  G4TwoVector v2(2601.2*mm,-417.4*mm);
  G4TwoVector v1(1777.6*mm,-433.5*mm);
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
  G4VSolid* steel_plate =  new G4ExtrudedSolid("SteelPlate",
					       vertexes,
					       size_z  / 2.0,
					       zero, 1.0,
					       zero, 1.0);

  //  DisplayVolume(steel_plate, hcalenvelope);
  volume_steel = steel_plate->GetCubicVolume()*n_scinti_plates;
  return steel_plate;
}

// Construct the envelope and the call the
// actual inner hcal construction
void
PHG4Prototype2OuterHcalDetector::Construct( G4LogicalVolume* logicWorld )
{
  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4VSolid* hcal_envelope_cylinder = new G4Tubs("OuterHcal_envelope_solid",  envelope_inner_radius-0.5*cm, envelope_outer_radius, envelope_z / 2., tan(-459.8/1770.9), fabs(tan(-459.8/1770.9)) +  12.4/180.* M_PI);
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
  ConstructOuterHcal(hcal_envelope_log);
  //  AddGeometryNode();
  return;
}

int
PHG4Prototype2OuterHcalDetector::ConstructOuterHcal(G4LogicalVolume* hcalenvelope)
{
  G4VSolid *steel_plate  = ConstructSteelPlate(hcalenvelope);
  G4LogicalVolume *steel_logical = new G4LogicalVolume(steel_plate, G4Material::GetMaterial(params->get_string_param("material")), "HcalInnerSteelPlate", 0, 0, 0);
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(true);
  visattchk->SetColour(G4Colour::Grey());
  steel_logical->SetVisAttributes(visattchk);
  //  G4AssemblyVolume *scinti_mother_logical = ConstructHcalScintillatorAssembly(hcalenvelope);
  double phi = 0.;
  double deltaphi = 2 * M_PI / n_scinti_plates;
  deltaphi = 2 * M_PI / 320.;
  ostringstream name;
  double middlerad = outer_radius - (outer_radius - inner_radius) / 2.;
  double shiftslat = fabs(scinti_tile_x_lower - scinti_tile_x_upper)/2.;
  for (int i = 0; i < n_scinti_plates; i++)
    //  for (int i = 0; i < 2; i++)
    {
      G4RotationMatrix *Rot = new G4RotationMatrix();
      double ypos = sin(phi) * middlerad;
      double xpos = cos(phi) * middlerad;
      // the center of the scintillator is not the center of the inner hcal
      // but depends on the tilt angle. Therefore we need to shift
      // the center from the mid point
      ypos += sin((-tilt_angle)/rad - phi)*shiftslat;
      xpos -= cos((-tilt_angle)/rad - phi)*shiftslat;
      Rot->rotateZ(phi * rad + tilt_angle);
      G4ThreeVector g4vec(xpos, ypos, 0);
      //      scinti_mother_logical->MakeImprint(hcalenvelope, g4vec, Rot, i, overlapcheck);
      Rot = new G4RotationMatrix();
      Rot->rotateZ(-phi * rad);
      name.str("");
      name << "OuterHcalSteel_" << i;
      steel_absorber_vec.insert(new G4PVPlacement(Rot, G4ThreeVector(0, 0, 0), steel_logical, name.str().c_str(), hcalenvelope, 0, i, overlapcheck));
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
  G4VSolid *bigtile = ConstructScintillatorBox(hcalenvelope);
  // eta->theta
  G4double delta_eta = params->get_double_param("scinti_eta_coverage") / n_scinti_tiles;
  G4double eta = 0;
  G4double theta;
  G4double x[4];
  G4double z[4];
  ostringstream name;
  double overhang = (scinti_tile_x - (outer_radius - inner_radius)) / 2.;
  double offset = 1 * cm + overhang; // add 1cm to make sure the G4ExtrudedSolid
  // is larger than the tile so we do not have
  // funny edge effects when overlapping vols
  double scinti_gap_neighbor = params->get_double_param("scinti_gap_neighbor")*cm;
  for (int i = 0; i < n_scinti_tiles; i++)
    {
      theta = M_PI / 2 - PHG4Utils::get_theta(eta); // theta = 90 for eta=0
      x[0] = inner_radius - overhang;
      z[0] = tan(theta) * inner_radius;
      x[1] = outer_radius + overhang; // since the tile is tilted, x is not at the outer radius but beyond
      z[1] = tan(theta) * outer_radius;
      eta += delta_eta;
      theta = M_PI / 2 - PHG4Utils::get_theta(eta); // theta = 90 for eta=0
      x[2] = inner_radius - overhang;
      z[2] =  tan(theta) * inner_radius;
      x[3] =  outer_radius + overhang; // since the tile is tilted, x is not at the outer radius but beyond
      z[3] = tan(theta) * outer_radius;
      // apply gap between scintillators
      z[0] += scinti_gap_neighbor / 2.;
      z[1] += scinti_gap_neighbor / 2.;
      z[2] -= scinti_gap_neighbor / 2.;
      z[3] -= scinti_gap_neighbor / 2.;
      Point_2 leftsidelow(z[0], x[0]);
      Point_2 leftsidehigh(z[1], x[1]);
      x[0] = inner_radius - offset;
      z[0] = x_at_y(leftsidelow, leftsidehigh, x[0]);
      x[1] = outer_radius + offset;
      z[1] = x_at_y(leftsidelow, leftsidehigh, x[1]);
      Point_2 rightsidelow(z[2], x[2]);
      Point_2 rightsidehigh(z[3], x[3]);
      x[2] = outer_radius + offset;
      z[2] = x_at_y(rightsidelow, rightsidehigh, x[2]);
      x[3] = inner_radius - offset;
      z[3] = x_at_y(rightsidelow, rightsidehigh, x[3]);


      vector<G4TwoVector> vertexes;
      for (int j = 0; j < 4; j++)
	{
	  G4TwoVector v(x[j], z[j]);
	  vertexes.push_back(v);
	}
      G4TwoVector zero(0, 0);

      G4VSolid *scinti =  new G4ExtrudedSolid("ScintillatorTile",
					      vertexes,
					      scinti_tile_thickness + 0.2 * mm,
					      zero, 1.0,
					      zero, 1.0);
      G4RotationMatrix *rotm = new G4RotationMatrix();
      rotm->rotateX(-90 * deg);
      name.str("");
      name << "scintillator_" << i << "_left";
      G4VSolid *scinti_tile =  new G4IntersectionSolid(name.str(), bigtile, scinti, rotm, G4ThreeVector(-(inner_radius + outer_radius) / 2., 0, 0));
      scinti_tiles_vec[i + n_scinti_tiles] = scinti_tile;
      rotm = new G4RotationMatrix();
      rotm->rotateX(90 * deg);
      name.str("");
      name << "scintillator_" << i << "_right";
      scinti_tile =  new G4IntersectionSolid(name.str(), bigtile, scinti, rotm, G4ThreeVector(-(inner_radius + outer_radius) / 2., 0, 0));
      scinti_tiles_vec[n_scinti_tiles - i - 1] =  scinti_tile;
    }
  // for (unsigned int i=0; i<scinti_tiles_vec.size(); i++)
  //   {
  //     if (scinti_tiles_vec[i])
  // 	 {
  // 	   DisplayVolume(scinti_tiles_vec[i],hcalenvelope );
  // 	 }
  //   }

  return;
}

double
PHG4Prototype2OuterHcalDetector::x_at_y(Point_2 &p0, Point_2 &p1, double yin)
{
  double xret = NAN;
  double x[2];
  double y[2];
  x[0] = CGAL::to_double(p0.x());
  y[0] = CGAL::to_double(p0.y());
  x[1] = CGAL::to_double(p1.x());
  y[1] = CGAL::to_double(p1.y());
  Line_2 l(p0, p1);
  double newx = fabs(x[0]) + fabs(x[1]);
  Point_2 p0new(-newx, yin);
  Point_2 p1new(newx, yin);
  Segment_2 s(p0new, p1new);
  CGAL::Object result = CGAL::intersection(l, s);
  if ( const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
    {
      xret = CGAL::to_double(ipoint->x());
    }
  else
    {
      cout << PHWHERE << " failed for y = " << y << endl;
      cout << "p0(x): " << CGAL::to_double(p0.x()) << ", p0(y): " <<  CGAL::to_double(p0.y()) << endl;
      cout << "p1(x): " << CGAL::to_double(p1.x()) << ", p1(y): " <<  CGAL::to_double(p1.y()) << endl;
      exit(1);
    }
  return xret;
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

