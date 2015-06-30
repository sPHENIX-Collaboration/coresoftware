#include "PHG4InnerHcalDetector.h"
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

static double no_overlap = 0.000;//15 * cm; // added safety margin against overlaps by using same boundary between volumes
PHG4InnerHcalDetector::PHG4InnerHcalDetector( PHCompositeNode *Node, const std::string &dnam, const int lyr  ):
  PHG4Detector(Node, dnam),
  inner_radius(1160 * mm),
  outer_radius(1360 * mm),
  size_z(1759.4 * 2 * mm),
  scinti_gap(8.5 * mm),
  tilt_angle(32.5 * deg),
  n_scinti_plates(5 * 64),
  scinti_tile_y(7 * mm),
  scinti_tile_z(size_z),
  envelope_inner_radius(1160 * mm - no_overlap),
  envelope_outer_radius(1360 * mm + no_overlap),
  envelope_z(size_z + no_overlap),
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
  scintilogicnameprefix("HcalInnerScinti")
{
}

PHG4InnerHcalDetector::~PHG4InnerHcalDetector()
{}

double
PHG4InnerHcalDetector::CalculateSteelAngularCoverage()
{
  //  double scinti_gap_projection = scinti_gap * cos(tilt_angle / rad);
  double scinti_gap_projection = (scinti_gap+2*mm)/cos(tilt_angle / rad);
   double cosphi = (2 * envelope_inner_radius * envelope_inner_radius - scinti_gap_projection * scinti_gap_projection) / (2 * envelope_inner_radius * envelope_inner_radius);
   double single_scinti_angcov = acos(cosphi);
  // double tanphi = scinti_gap_projection/envelope_inner_radius;
  // double single_scinti_angcov = atan(tanphi);

  double total_scinti_coverage = single_scinti_angcov * n_scinti_plates;
  double total_steel_coverage = 2 * M_PI - total_scinti_coverage;
  double single_steel_coverage = total_steel_coverage / n_scinti_plates;
  return single_steel_coverage;
}

//_______________________________________________________________
//_______________________________________________________________
int
PHG4InnerHcalDetector::IsInInnerHcal(G4VPhysicalVolume * volume) const
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
  // if (absorberactive)
  //   {
  //     if (steel_absorber_vec.find(volume) != steel_absorber_vec.end())
  // 	{
  // 	  return -1;
  // 	}
  //   }
  // if (volume->GetName().find(scintilogicnameprefix) != string::npos)
  //   {
  //     return 1;
  //   }
  return 0;
}

G4VSolid*
PHG4InnerHcalDetector::ConstructScintillatorBox(G4LogicalVolume* hcalenvelope)
{
  double mid_radius = inner_radius + (outer_radius-inner_radius)/2.;
  Point_2 p_in_1(mid_radius,0); // center of scintillator
  double angle_mid_scinti = M_PI/2. - fabs(tilt_angle);
  cout << "scinti center: x " << CGAL::to_double(p_in_1.x()) << ", y: " <<  CGAL::to_double(p_in_1.y()) << endl;
  // x coordinate of end of center vertical
  double xcoord = scinti_tile_y/2. * cos(angle_mid_scinti*rad) + mid_radius;
  double ycoord =   scinti_tile_y/2. * sin(angle_mid_scinti*rad) + 0;
  Point_2 p_upperedge(xcoord,ycoord);
  Line_2 s2(p_in_1,p_upperedge); // center vertical

  Line_2 perp =  s2.perpendicular(p_upperedge);
  Point_2 sc1(outer_radius,0), sc2(0,outer_radius),sc3(-outer_radius,0);
  Circle_2 outer_circle(sc1,sc2,sc3);
  vector< CGAL::Object > res;
  CGAL::intersection(outer_circle, perp, std::back_inserter(res));
  Point_2 upperright;
  vector< CGAL::Object >::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); iter++)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) >  CGAL::to_double(p_upperedge.x()))
	    {
	      cout << "std::pair<Circular_arc_point_2<Circular_k>, unsigned>" << endl;
	      //	      cout << "intersect: " << point->first << ", n: " << point->second << endl;
	      cout << "upper right x: " << CGAL::to_double(point->first.x()) << ", y: " << CGAL::to_double(point->first.y()) << endl;
	      double deltax = CGAL::to_double(point->first.x())-CGAL::to_double(p_upperedge.x());
	      double deltay = CGAL::to_double(point->first.y())-CGAL::to_double(p_upperedge.y());
	      cout << "dist: " << sqrt(deltax*deltax+deltay*deltay) << endl;
	      // sqrt(deltax*deltax+deltay*deltay) is distance from scintilator center to outer edge
	      // the scintillator is twice as long
	      scinti_tile_x = 2*sqrt(deltax*deltax+deltay*deltay); // 
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      upperright = pntmp;
	    }
	}
      else
	{
	  cout << "CGAL::Object type not pair..." << endl;
	}
    }
  cout << "scinti_tile_x : " << scinti_tile_x << ", scinti_tile_y: " << scinti_tile_y
       << ", scinti_tile_z: " << scinti_tile_z << endl;
  G4VSolid* scintibox =  new G4Box("ScintiTile", scinti_tile_x / 2., scinti_tile_y / 2., scinti_tile_z / 2.);
 
  return scintibox;
}

G4VSolid*
PHG4InnerHcalDetector::ConstructScintillatorBoxA(G4LogicalVolume* hcalenvelope)
{
  // procedure:
  // we need to calculate the dimension of the box which is tilted inside the cylinder
  // we start with the inner envelope at x=0, create a line with
  // the tilt angle, create a line perpendicular to it at the intersection, get the point
  // on that line at the thickness of the scintillator (upper left cornet), then create a line 
  // with the tilt angle at it and find the intersection with the outer envelope, that gives us the
  // length of the scintillator
  Point_2 p_in_1(inner_radius,0); // point where the left scintilator corner touches the envelope
  Point_2 py(inner_radius+1,tan(tilt_angle*rad)); // random point on scintilator edge to get the line
  Line_2 s2(p_in_1,py); // edge of the scintillator
  Line_2 perp =  s2.perpendicular(p_in_1); // perpendicular line the left scintilator edge
  Point_2 sc1(inner_radius+scinti_tile_y,0), sc2(inner_radius-scinti_tile_y,0),sc3(inner_radius,scinti_tile_y);
  Circle_2 scithk(sc1,sc2,sc3);
  vector< CGAL::Object > res;
  CGAL::intersection(scithk, perp, std::back_inserter(res));
  Point_2 upperleft;
  vector< CGAL::Object >::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); ++iter)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
   	{
  	  if (CGAL::to_double(point->first.x()) > CGAL::to_double(p_in_1.x()))
	    {
	      cout << "std::pair<Circular_arc_point_2<Circular_k>, unsigned>" << endl;
	      cout << "upper left x: " << CGAL::to_double(point->first.x()) << ", y: " << CGAL::to_double(point->first.y()) << endl;
	      double deltax = CGAL::to_double(point->first.x())-CGAL::to_double(p_in_1.x());
	      double deltay = CGAL::to_double(point->first.y())-CGAL::to_double(p_in_1.y());
	      cout << "dist: " << sqrt(deltax*deltax+deltay*deltay) << endl;
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      upperleft = pntmp;
	    }
	}
      else
	{
	  cout << "CGAL::Object type not pair..." << endl;
	}
    }
  Line_2 upperline  = perp.perpendicular(upperleft);
  Point_2 outrad1(outer_radius,0), outrad2(0,outer_radius), outrad3(-outer_radius,0);
  Circle_2 outerradius(outrad1,outrad2,outrad3);
  res.clear();
  CGAL::intersection(outerradius, upperline, std::back_inserter(res));
  Point_2 upperright;
  for (iter = res.begin(); iter != res.end(); iter++)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
  	{
  	  if (CGAL::to_double(point->first.x()) > CGAL::to_double(upperleft.x()))
  	    {
  	      cout << "std::pair<Circular_arc_point_2<Circular_k>, unsigned>" << endl;
  	      //	      cout << "intersect: " << point->first << ", n: " << point->second << endl;
  	      cout << "upper right x: " << CGAL::to_double(point->first.x()) << ", y: " << CGAL::to_double(point->first.y()) << endl;
  	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
  	      upperright = pntmp;
  	      double deltax = CGAL::to_double(point->first.x())-CGAL::to_double(upperleft.x());
  	      double deltay = CGAL::to_double(point->first.y())-CGAL::to_double(upperleft.y());
  	      cout << "scintilen: " << sqrt(deltax*deltax+deltay*deltay) << endl;
  	      scinti_tile_x = sqrt(deltax*deltax+deltay*deltay);
  	    }
  	}
      else
  	{
  	  cout << "CGAL::Object type not pair..." << endl;
  	}
    }

  cout << "scinti_tile_x : " << scinti_tile_x << ", scinti_tile_y: " << scinti_tile_y
       << ", scinti_tile_z: " << scinti_tile_z << endl;
  G4VSolid* scintibox =  new G4Box("ScintiTile", scinti_tile_x / 2., scinti_tile_y / 2., scinti_tile_z / 2.);
 
  return scintibox;
     }


G4VSolid*
PHG4InnerHcalDetector::ConstructSteelPlate(G4LogicalVolume* hcalenvelope)
{
  //   A                  C
  //   *------------------*
  //    \                 |
  //     *M               |
  //      \               |
  //       *--------------*
  //       B              D
  // procedure:
  // find the intersection of the inner radius circle with the radius at 1/2 the
  // coverage angle (M: where the steel plate will touch the inner envelope)
  // construct the tangent at that point (line perpendicular to radius)
  // and get intersections with radii at 0 (A) and coverage angle (B) which are
  // the left corners of steel plate
  // construct lines from corners with tilt angle and get intersections of those
  // with outer radius circle to find right corners of steel plate (C,D)
  Point_2 pnull(0, 0);
  Point_2 p_in_1(0, inner_radius), p_in_2(inner_radius, 0), p_in_3(-inner_radius, 0);
  Circle_2 c_inner(p_in_1, p_in_2, p_in_3);
  Point_2 p_out_1(0, outer_radius), p_out_2(outer_radius, 0), p_out_3(-outer_radius, 0);
  Circle_2 c_outer(p_out_1, p_out_2, p_out_3);
  // find M
  Point_2 py(-1, tan(single_steel_angular_coverage / 2.));
  Line_2 s2(pnull, py);
  vector< CGAL::Object > res;
  CGAL::intersection(c_inner, s2, std::back_inserter(res));
  Point_2 intersection_center_front;
  vector< CGAL::Object >::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); iter++)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) > 0)
	    {
	      cout << "std::pair<Circular_arc_point_2<Circular_k>, unsigned>" << endl;
	      cout << "intersect: " << point->first << ", n: " << point->second << endl;
	      cout << "tangente x: " << CGAL::to_double(point->first.x()) << ", tangente y: " << CGAL::to_double(point->first.y()) << endl;
	    }
	  Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	  intersection_center_front = pntmp;
	}
      else
	{
	  cout << "CGAL::Object type not pair..." << endl;
	}
    }
  CGAL::Object result;
  Line_2 tangent = s2.perpendicular(intersection_center_front); // tangent at intersection
  // CGAL only knows intersections of lines (infinite) with segments (lines with begin/end point)
  // so construct a segment long enough to intersect with our tangent, starting at 0/0
  Point_2 pnullfar(2 * inner_radius, 0);
  Segment_2 segnull(pnull, pnullfar);
  Point_2 upperleftcorner;
  result = CGAL::intersection(segnull, tangent);
  if (CGAL::assign(upperleftcorner, result))
    {
      cout << "upper left x: " << CGAL::to_double(upperleftcorner.x())
	   << ", upper left y: " << CGAL::to_double(upperleftcorner.y())
	   << endl;
    }
  else
    {
      cout << "no usable intersection for upper left corner" << endl;
      cout << "pnull: x " <<  CGAL::to_double(pnull.x())
	   << ", y " << CGAL::to_double(pnull.y())
	   << endl;
      cout << "pnullfar: x " <<  CGAL::to_double(pnullfar.x())
	   << ", y " << CGAL::to_double(pnullfar.y())
	   << endl;
      exit(1);
    }
  Point_2 plowcornerfar(2 * inner_radius, -(2 * inner_radius * tan(single_steel_angular_coverage)));
  Segment_2 seglowcorner(pnull, plowcornerfar);
  Point_2 lowerleftcorner;
  result = CGAL::intersection(seglowcorner, tangent);
  if (CGAL::assign(lowerleftcorner, result))
    {
      cout << "lower left x : " << CGAL::to_double(lowerleftcorner.x())
	   << ", lower left y: " << CGAL::to_double(lowerleftcorner.y())
	   << endl;
      // handle the point intersection case.

    }
  else
    {
      cout << "no usable intersection for lower left cornet" << endl;
      exit(1);
    }
  // now for the right corners, construct line from left corners with
  // slope of tilt angle and intersect with circle with outer radius
  // upper left cornet first
  //			    CGAL::to_double(upperleftcorner.y()) - (abs(tilt_angle) / tilt_angle)*tan(tilt_angle));
  Point_2 upperleftcorner_1(CGAL::to_double(upperleftcorner.x()) + 1,
			    CGAL::to_double(upperleftcorner.y()) - (boost::math::sign(tilt_angle)*tan(tilt_angle)));
  Line_2 upperborder(upperleftcorner, upperleftcorner_1);
  res.clear(); // just clear the content from the last intersection search
  CGAL::intersection(c_outer, upperborder, std::back_inserter(res));
  Point_2 upperrightcorner;
  for (iter = res.begin(); iter != res.end(); iter++)
    {
      CGAL::Object obj = *iter;

      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) > 0)
	    {
	      cout << "upperright x: " << CGAL::to_double(point->first.x()) << ", upperright y: " << CGAL::to_double(point->first.y()) << endl;
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      upperrightcorner = pntmp;
	    }
	}
    }
  // lower right corner
  Point_2 lowerleftcorner_1(CGAL::to_double(lowerleftcorner.x()) + 1,
			    CGAL::to_double(lowerleftcorner.y()) - (boost::math::sign(tilt_angle)*tan(tilt_angle+boost::math::sign(tilt_angle)*2*M_PI/n_scinti_plates)));
  Line_2 lowerborder(lowerleftcorner, lowerleftcorner_1 );
  Point_2 lowerrightcorner;
  res.clear();
  CGAL::intersection(c_outer, lowerborder, std::back_inserter(res));
  for (iter = res.begin(); iter != res.end(); iter++)
    {
      CGAL::Object obj = *iter;

      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) > 0)
	    {
	      cout << "lowerright x: " << CGAL::to_double(point->first.x()) << ", lowerright y: " << CGAL::to_double(point->first.y()) << endl;
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      lowerrightcorner = pntmp;
	    }
	}
    }
  G4TwoVector v1(CGAL::to_double(upperleftcorner.x()), CGAL::to_double(upperleftcorner.y()));
  G4TwoVector v2(CGAL::to_double(upperrightcorner.x()), CGAL::to_double(upperrightcorner.y()));
  G4TwoVector v3(CGAL::to_double(lowerrightcorner.x()), CGAL::to_double(lowerrightcorner.y()));
  G4TwoVector v4(CGAL::to_double(lowerleftcorner.x()), CGAL::to_double(lowerleftcorner.y()));
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
  return steel_plate;
}

void
PHG4InnerHcalDetector::ConstructScintillator(G4LogicalVolume *hcalenvelope)
{
  scinti_tile_x = outer_radius - inner_radius;
  G4VSolid* scinti_tile =  new G4Box("ScintiTile", scinti_tile_x / 2., scinti_tile_y / 2., scinti_tile_z / 2.);
  DisplayVolume(scinti_tile, hcalenvelope);
  return;
}


// 
void
PHG4InnerHcalDetector::Construct( G4LogicalVolume* logicWorld )
{
  // first calculate some parameters which depend on
  // macro-changeable configurations
  single_steel_angular_coverage = CalculateSteelAngularCoverage();
  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4VSolid* hcal_envelope_cylinder = new G4Tubs("InnerHcal_envelope_solid",  envelope_inner_radius, envelope_outer_radius, envelope_z / 2., 0, 2 * M_PI);
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
  new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(place_in_x, place_in_y, place_in_z)), hcal_envelope_log, "InnerHcalEnvelope", logicWorld, 0, false, overlapcheck);
  ConstructInnerHcal(hcal_envelope_log);
  //  AddGeometryNode();
  return;
}

int
PHG4InnerHcalDetector::ConstructInnerHcal(G4LogicalVolume* hcalenvelope)
{
  G4VSolid *steel_plate  = ConstructSteelPlate(hcalenvelope);
  G4LogicalVolume *steel_logical = new G4LogicalVolume(steel_plate, G4Material::GetMaterial("G4_Fe"), "HcalOuterSteelPlate", 0, 0, 0);
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  visattchk->SetColour(G4Colour::Cyan());
  steel_logical->SetVisAttributes(visattchk);
  G4VSolid *scintibox  = ConstructScintillatorBox(hcalenvelope);
  G4LogicalVolume *scintibox_logical = new G4LogicalVolume(scintibox, G4Material::GetMaterial("G4_POLYSTYRENE"), "HcalScintiBox", 0, 0, 0);
  visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  visattchk->SetColour(G4Colour::Red());
  scintibox_logical->SetVisAttributes(visattchk);
  double thickness = single_steel_angular_coverage;
  thickness += scinti_gap;
  double phi = 0;
  //  double phioff = (2*M_PI/n_scinti_plates)*2+0.005;
  //double deltaphi = acos((2 * envelope_inner_radius * envelope_inner_radius - thickness * thickness) / (2 * envelope_inner_radius * envelope_inner_radius));
  double deltaphi = 2*M_PI/n_scinti_plates;
  double xpos = 0;
  double ypos = 0;
  ostringstream name;
  double middlerad = outer_radius - (outer_radius-inner_radius)/2.;
  for (int i = 0; i < n_scinti_plates; i++)
    {
      G4RotationMatrix *Rot = new G4RotationMatrix();
      Rot->rotateZ(-tilt_angle * rad-phi * rad);
      ypos = sin(phi)*middlerad;
      xpos = cos(phi)*middlerad;
      name.str("");
      name << "ScintiBox_" << i;
      new G4PVPlacement(Rot, G4ThreeVector(xpos, ypos, 0), scintibox_logical, name.str().c_str(), hcalenvelope, 0, i, overlapcheck);
      // Rot = new G4RotationMatrix();
      // Rot->rotateZ(-phi * rad);
      // name.str("");
      // name << "HcalSteel_" << i;
      //      new G4PVPlacement(Rot, G4ThreeVector(0, 0, 0), steel_logical, name.str().c_str(), hcalenvelope, 0, i, overlapcheck);
      phi += deltaphi;
    }
  //  ConstructScintillator(hcalenvelope);
  return 0;
}


int
PHG4InnerHcalDetector::DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix *rotm )
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
PHG4InnerHcalDetector::AddGeometryNode()
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
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeomv3(inner_radius / cm, (place_in_z - size_z / 2.) / cm, (place_in_z + size_z / 2.) / cm, (outer_radius - inner_radius) / cm, n_scinti_plates,  tilt_angle / rad, 0);
      geo->AddLayerGeom(layer, mygeom);
      geo->identify();
    }
}
