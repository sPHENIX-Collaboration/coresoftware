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

//static double no_overlap = 0.00015 * cm; // added safety margin against overlaps by using same boundary between volumes
// face touches the boundary instead of the corner, subtracting 1 permille from the total
// scintilator length takes care of this
static double subtract_from_scinti_x = 0.1*mm;

PHG4OuterHcalDetector::PHG4OuterHcalDetector( PHCompositeNode *Node, PHG4OuterHcalParameters *parameters,const std::string &dnam):
  PHG4Detector(Node, dnam),
  params(parameters),
  envelope_inner_radius(params->inner_radius),
  envelope_outer_radius(params->outer_radius),
  envelope_z(params->size_z),
  scinti_tile_x(NAN),
  scinti_tile_x_lower(NAN),
  scinti_tile_x_upper(NAN),
  scinti_tile_z(params->size_z),
  steel_rectangle_plate_x(657.2*mm),
  steel_plate_x(828.7*mm),
  steel_plate_z(3049.1*2*mm),
  n_steel_plates(320),
  scinti_tile_x_old(821.1*mm),
  scinti_tile_z_old(steel_plate_z),
  scinti_eta_coverage(1.1),
  etacutline(0.8),
  cutbox_x((steel_plate_x - steel_rectangle_plate_x*mm)*2),// twice the size we need to cut to make geo easier
  // trapezoid twice the size of what we need
  cuttrapezoid_x(cutbox_x),
  layer(0),
  scintilogicnameprefix("HcalOuterScinti"),
  field_setup(NULL)
{
  scinti_tiles_vec.assign(2*params->n_scinti_tiles,static_cast<G4VSolid *>(NULL));

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
  G4double steel_surface = total_surface - params->scinti_gap;
  // 1*mm/320. is a fudge factor, if using the calculated numbers
  // I get a 130um overlap when putting the last scintilator in
  // We have an air gap of (8.5-7)/2mm = 0.75mm + 0.13mm = 0.88mm ~ 1mm too much
  // We divide this by the number of panels to adjust the thickness of each steel plate
  G4double fudge_subtract = 1.*mm/n_steel_plates;
  steel_plate_yin = steel_surface/cos(params->tilt_angle/rad)-fudge_subtract;
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
  steel_surface = total_surface - params->scinti_gap;
  steel_plate_yout = steel_surface/cos(params->tilt_angle/rad)-fudge_subtract;

  // allocate memory for scintillator plates
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
  if (params->absorberactive)
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

G4VSolid*
PHG4OuterHcalDetector::ConstructScintillatorBox(G4LogicalVolume* hcalenvelope)
{
  double mid_radius = params->inner_radius + (params->outer_radius - params->inner_radius) / 2.;
  Point_2 p_in_1(mid_radius, 0); // center of scintillator

  // length of upper edge (middle till outer circle intersect
  // x/y coordinate of end of center vertical
  double xcoord  = params->scinti_tile_thickness / 2. *sin(fabs(params->tilt_angle)/rad) +  mid_radius;
  double ycoord  =   params->scinti_tile_thickness / 2. * cos(fabs(params->tilt_angle)/rad) + 0;
  Point_2 p_upperedge(xcoord, ycoord);
  Line_2 s2(p_in_1, p_upperedge); // center vertical

  Line_2 perp =  s2.perpendicular(p_upperedge);
  Point_2 sc1(params->outer_radius, 0), sc2(0, params->outer_radius), sc3(-params->outer_radius, 0);
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
  xcoord = mid_radius - params->scinti_tile_thickness / 2. * sin(fabs(params->tilt_angle)/rad);
  ycoord  = 0 -  params->scinti_tile_thickness / 2. * cos(fabs(params->tilt_angle)/rad);
  Point_2 p_loweredge(xcoord, ycoord);
  Line_2 s3(p_in_1, p_loweredge);
  Line_2 l_lower = s3.perpendicular(p_loweredge);
  Point_2 ic1(params->inner_radius, 0), ic2(0, params->inner_radius), ic3(-params->inner_radius, 0);
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
  G4VSolid* scintibox =  new G4Box("ScintiTile", scinti_tile_x / 2., params->scinti_tile_thickness / 2., scinti_tile_z / 2.);

  return scintibox;
}

G4VSolid*
PHG4OuterHcalDetector::ConstructSteelPlate(G4LogicalVolume* hcalenvelope)
{
  // calculate steel plate on top of the scinti box. Lower edge is the upper edge of
  // the scintibox + 1/2 the airgap
  double mid_radius = params->inner_radius + (params->outer_radius - params->inner_radius) / 2.;
  // first the lower edge, just like the scinti box, just add the air gap
  // and calculate intersection of edge with inner and outer radius.
  Point_2 p_in_1(mid_radius, 0); // center of lower scintillator
  double angle_mid_scinti = M_PI / 2. + params->tilt_angle / rad;
  double xcoord = params->scinti_gap / 2. * cos(angle_mid_scinti / rad) + mid_radius;
  double ycoord =   params->scinti_gap / 2. * sin(angle_mid_scinti / rad) + 0;
  Point_2 p_loweredge(xcoord, ycoord);
  Line_2 s2(p_in_1, p_loweredge); // center vertical
  Line_2 perp =  s2.perpendicular(p_loweredge); // that is the lower edge of the steel plate
  Point_2 sc1(params->inner_radius, 0), sc2(0, params->inner_radius), sc3(-params->inner_radius, 0);
  Circle_2 inner_circle(sc1, sc2, sc3);
  vector< CGAL::Object > res;
  CGAL::intersection(inner_circle, perp, std::back_inserter(res));
  Point_2 lowerleft;
  vector< CGAL::Object >::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); ++iter)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) > 0)
	    {
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      lowerleft = pntmp;
	    }
	}
      else
	{
	  cout << "CGAL::Object type not pair..." << endl;
	}
    }
  Point_2 so1(params->outer_radius, 0), so2(0, params->outer_radius), so3(-params->outer_radius, 0);
  Circle_2 outer_circle(so1, so2, so3);
  res.clear(); // just clear the content from the last intersection search
  CGAL::intersection(outer_circle, perp, std::back_inserter(res));
  Point_2 lowerright;
  for (iter = res.begin(); iter != res.end(); ++iter)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) >  CGAL::to_double(p_loweredge.x()))
	    {
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      lowerright = pntmp;
	    }
	}
      else
	{
	  cout << "CGAL::Object type not pair..." << endl;
	}
    }
  // now we have the lower left and rigth corner, now find the upper edge
  // find the center of the upper scintilator


  double phi_midpoint = 2 * M_PI / params->n_scinti_plates;
  double xmidpoint = cos(phi_midpoint) * mid_radius;
  double ymidpoint = sin(phi_midpoint) * mid_radius;
  // angle of perp line at center of scintillator
  angle_mid_scinti = (M_PI / 2. - phi_midpoint) - (M_PI / 2. + params->tilt_angle / rad);
  double xcoordup = xmidpoint - params->scinti_gap / 2. * sin(angle_mid_scinti / rad);
  double ycoordup = ymidpoint - params->scinti_gap / 2. * cos(angle_mid_scinti / rad);
  Point_2 upperleft;
  Point_2 upperright;
  Point_2 mid_upperscint(xmidpoint, ymidpoint);
  Point_2 p_upperedge(xcoordup, ycoordup);
  {
    Line_2 sup(mid_upperscint, p_upperedge); // center vertical
    Line_2 perp =  sup.perpendicular(p_upperedge); // that is the upper edge of the steel plate
    Point_2 sc1(params->inner_radius, 0), sc2(0, params->inner_radius), sc3(-params->inner_radius, 0);
    Circle_2 inner_circle(sc1, sc2, sc3);
    vector< CGAL::Object > res;
    CGAL::intersection(inner_circle, perp, std::back_inserter(res));
    vector< CGAL::Object >::const_iterator iter;
    double pxmax = 0.;
    for (iter = res.begin(); iter != res.end(); ++iter)
      {
	CGAL::Object obj = *iter;
	if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	  {
	    if (CGAL::to_double(point->first.x()) > pxmax)
	      {
		pxmax = CGAL::to_double(point->first.x());
		Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
		upperleft = pntmp;
	      }
	  }
	else
	  {
	    cout << "CGAL::Object type not pair..." << endl;
	  }
      }
    Point_2 so1(params->outer_radius, 0), so2(0, params->outer_radius), so3(-params->outer_radius, 0);
    Circle_2 outer_circle(so1, so2, so3);
    res.clear(); // just clear the content from the last intersection search
    CGAL::intersection(outer_circle, perp, std::back_inserter(res));
    for (iter = res.begin(); iter != res.end(); ++iter)
      {
	CGAL::Object obj = *iter;
	if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	  {
	    if (CGAL::to_double(point->first.x()) >  CGAL::to_double(p_loweredge.x()))
	      {
		Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
		upperright = pntmp;
	      }
	  }
	else
	  {
	    cout << "CGAL::Object type not pair..." << endl;
	  }
      }
  }
  // the left corners are on a secant with the inner boundary, they need to be shifted
  // to be a tangent at the center
  ShiftSecantToTangent(lowerleft, upperleft, upperright, lowerright);
  G4TwoVector v1(CGAL::to_double(upperleft.x()), CGAL::to_double(upperleft.y()));
  G4TwoVector v2(CGAL::to_double(upperright.x()), CGAL::to_double(upperright.y()));
  G4TwoVector v3(CGAL::to_double(lowerright.x()), CGAL::to_double(lowerright.y()));
  G4TwoVector v4(CGAL::to_double(lowerleft.x()), CGAL::to_double(lowerleft.y()));
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(v1);
  vertexes.push_back(v2);
  vertexes.push_back(v3);
  vertexes.push_back(v4);
  G4TwoVector zero(0, 0);
  G4VSolid* steel_plate =  new G4ExtrudedSolid("SteelPlate",
					       vertexes,
					       params->size_z  / 2.0,
					       zero, 1.0,
					       zero, 1.0);

  //  DisplayVolume(steel_plate, hcalenvelope);
 
 return steel_plate;
}

void
PHG4OuterHcalDetector::ShiftSecantToTangent(Point_2 &lowleft, Point_2 &upleft, Point_2 &upright, Point_2 &lowright)
{
  Line_2 secant(lowleft, upleft);
  Segment_2 upedge(upleft, upright);
  Segment_2 lowedge(lowleft, lowright);
  double xmid = (CGAL::to_double(lowleft.x()) + CGAL::to_double(upleft.x())) / 2.;
  double ymid = (CGAL::to_double(lowleft.y()) + CGAL::to_double(upleft.y())) / 2.;
  Point_2 midpoint(xmid, ymid);
  Line_2 sekperp = secant.perpendicular(midpoint);
  Point_2 sc1(params->inner_radius, 0), sc2(0, params->inner_radius), sc3(-params->inner_radius, 0);
  Circle_2 inner_circle(sc1, sc2, sc3);
  vector< CGAL::Object > res;
  CGAL::intersection(inner_circle, sekperp, std::back_inserter(res));
  vector< CGAL::Object >::const_iterator iter;
  double pxmax = 0.;
  Point_2 tangtouch;
  for (iter = res.begin(); iter != res.end(); ++iter)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) > pxmax)
	    {
	      pxmax = CGAL::to_double(point->first.x());
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      tangtouch = pntmp;
	    }
	}
      else
	{
	  cout << "CGAL::Object type not pair..." << endl;
	}
    }
  Line_2 leftside = sekperp.perpendicular(tangtouch);
  CGAL::Object result = CGAL::intersection(upedge, leftside);
  if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
    {
      upleft = *ipoint;
    }
  result = CGAL::intersection(lowedge, leftside);
  if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
    {
      lowleft = *ipoint;
    }
  return;
}

void
PHG4OuterHcalDetector::Construct( G4LogicalVolume* logicWorld )
{
  field_setup = new PHG4OuterHcalFieldSetup(
  params->n_scinti_plates,/*G4int steelPlates*/
  params->scinti_gap, /*G4double scintiGap*/
  params->tilt_angle);/*G4double tiltAngle*/


  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4VSolid* hcal_envelope_cylinder = new G4Tubs("OuterHcal_envelope_solid",  envelope_inner_radius, envelope_outer_radius, envelope_z/2.,0,2*M_PI);
  G4LogicalVolume* hcal_envelope_log =  new G4LogicalVolume(hcal_envelope_cylinder, Air, G4String("OuterHcal_envelope"), 0, 0, 0);
  G4VisAttributes* hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Magenta());
  hcal_envelope_log->SetVisAttributes(hcalVisAtt);
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(params->x_rot);
  hcal_rotm.rotateY(params->y_rot);
  hcal_rotm.rotateZ(params->z_rot);
  new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(params->place_in_x, params->place_in_y, params->place_in_z)), hcal_envelope_log, "OuterHcal", logicWorld, 0, false, overlapcheck);
  ConstructOuterHcal(hcal_envelope_log);
  AddGeometryNode();
  return;
}

int
PHG4OuterHcalDetector::ConstructOuterHcal(G4LogicalVolume* hcalenvelope)
{
  ConsistencyCheck();
  SetTiltViaNcross(); // if number of crossings is set, use it to determine tilt
  CheckTiltAngle(); // die if the tilt angle is out of range
  G4VSolid *steel_plate =  ConstructSteelPlate(hcalenvelope);
  //   DisplayVolume(steel_plate_4 ,hcalenvelope);
  G4LogicalVolume *steel_logical = new G4LogicalVolume(steel_plate, G4Material::GetMaterial(params->material), "HcalOuterSteelPlate", 0, 0, 0);
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  visattchk->SetColour(G4Colour::Cyan());
  steel_logical->SetVisAttributes(visattchk);
  G4AssemblyVolume *scinti_mother_logical = ConstructHcalScintillatorAssembly(hcalenvelope);
  double phi = 0;
  double deltaphi = 2 * M_PI / params->n_scinti_plates;
  ostringstream name;
  double middlerad = params->outer_radius - (params->outer_radius - params->inner_radius) / 2.;
  double shiftslat = fabs(scinti_tile_x_lower - scinti_tile_x_upper)/2.;
  //  for (int i = 0; i < params->n_scinti_plates; i++)
  for (int i = 0; i < 1; i++)
    {
      G4RotationMatrix *Rot = new G4RotationMatrix();
      double ypos = sin(phi) * middlerad;
      double xpos = cos(phi) * middlerad;
      // the center of the scintillator is not the center of the inner hcal
      // but depends on the tilt angle. Therefore we need to shift
      // the center from the mid point
      ypos += sin((-params->tilt_angle)/rad - phi)*shiftslat;
      xpos -= cos((-params->tilt_angle)/rad - phi)*shiftslat;
      Rot->rotateZ(phi * rad + params->tilt_angle);
      G4ThreeVector g4vec(xpos, ypos, 0);
      scinti_mother_logical->MakeImprint(hcalenvelope, g4vec, Rot, i, overlapcheck);
      Rot = new G4RotationMatrix();
      Rot->rotateZ(-phi * rad);
      name.str("");
      name << "OuterHcalSteel_" << i;
      steel_absorber_vec.insert(new G4PVPlacement(Rot, G4ThreeVector(0, 0, 0), steel_logical, name.str().c_str(), hcalenvelope, 0, i, overlapcheck));
      phi += deltaphi;
    }
  hcalenvelope->SetFieldManager(
      field_setup -> get_Field_Manager_Gap(),
      false);

  steel_logical->SetFieldManager(
      field_setup -> get_Field_Manager_Iron(),
      true);
  return 0;
}

int
PHG4OuterHcalDetector::ConstructOuterHcal_A(G4LogicalVolume* hcalenvelope)
{
  ConsistencyCheck();
  SetTiltViaNcross(); // if number of crossings is set, use it to determine tilt
  CheckTiltAngle(); // die if the tilt angle is out of range
  G4VSolid *steel_plate =  ConstructSteelPlate(hcalenvelope);
  //   DisplayVolume(steel_plate_4 ,hcalenvelope);
  G4LogicalVolume *steel_logical = new G4LogicalVolume(steel_plate, G4Material::GetMaterial(params->material), "HcalOuterSteelPlate", 0, 0, 0);
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  visattchk->SetColour(G4Colour::Cyan());
  steel_logical->SetVisAttributes(visattchk);
  G4AssemblyVolume *scinti_mother_logical = ConstructHcalScintillatorAssembly(hcalenvelope);
  // visattchk = new G4VisAttributes();
  // visattchk->SetVisibility(true);
  // visattchk->SetForceSolid(false);
  // visattchk->SetColour(G4Colour::Red());
  // scinti_mother_logical->SetVisAttributes(visattchk);

  double thickness = steel_plate_yin * cos(params->tilt_angle / rad);
  thickness += params->scinti_gap;
  double deltaphi = acos((2 * envelope_inner_radius * envelope_inner_radius - thickness * thickness) / (2 * envelope_inner_radius * envelope_inner_radius));
  double phi = 0;
  ostringstream name;
  for (int i = 0; i < n_steel_plates; i++)
    {
      G4RotationMatrix *Rot = new G4RotationMatrix();
      name.str("");
      name << "HcalOuterScintiMother_" << i;
      Rot->rotateZ(phi * rad - params->tilt_angle);
      double xpos_scinti = (envelope_inner_radius+(scinti_tile_x_old)/2.) * cos(phi);
      double ypos_scinti = (envelope_inner_radius+(scinti_tile_x_old)/2.) * sin(phi);
      G4ThreeVector g4vec(xpos_scinti, ypos_scinti, 0);
      scinti_mother_logical->MakeImprint(hcalenvelope,g4vec,Rot,i,overlapcheck);
      Rot = new G4RotationMatrix();
      Rot->rotateZ(-phi * rad + params->tilt_angle);
      name.str("");
      name << "HcalOuterSteelPlate" << i;
      // start at the same position as the scintillator tiles. Naturally G4 has a different center for
      // rotating/translating the objects - the reference for the extruded solid seems to be
      // the upper left corner, the G4Box for the scinitllator has the center as reference
      // now shift this into the middle of the gap between the scintillator tiles at the inner radius
      // using the tilt angle and the angle of the slat
      double xpos = xpos_scinti;
      double ypos = ypos_scinti;
      // ypos += ((envelope_outer_radius-envelope_inner_radius)/2.)*sin(-phi/rad + params->tilt_angle/rad);
      // xpos -= ((envelope_outer_radius-envelope_inner_radius)/2.)*cos(-phi/rad + params->tilt_angle/rad);
      ypos += (scinti_tile_x_old/2.)*sin(-phi/rad + params->tilt_angle/rad);
      xpos -= (scinti_tile_x_old/2.)*cos(-phi/rad + params->tilt_angle/rad);
      // now shift the steel extruded shape down by half the width of the scintillator + gap
      // using the tilt angle and the angle of the slat
      xpos -= sin(-phi/rad + params->tilt_angle/rad)* (params->scinti_tile_thickness+(params->scinti_gap-params->scinti_tile_thickness)/2.)/2.;
      ypos -= cos(-phi/rad + params->tilt_angle/rad)* (params->scinti_tile_thickness+(params->scinti_gap-params->scinti_tile_thickness)/2.)/2.;
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
  G4double lowerright_x = steel_plate_x+sin(params->tilt_angle/rad)*steel_plate_yout;
  G4double lowerright_y = -cos(params->tilt_angle/rad)*steel_plate_yout;
  G4double lowerleft_x = sin(params->tilt_angle/rad)*steel_plate_yin;
  G4double lowerleft_y = -cos(params->tilt_angle/rad)*steel_plate_yin;
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
  // accomodate for the 12 deg params->tilt_angle)
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
  double yshift = tan(params->tilt_angle/rad)*(steel_plate_x - steel_rectangle_plate_x*mm);
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
  G4VSolid* scinti_tile =  new G4Box("ScintiTile",scinti_tile_x_old/2.,params->scinti_tile_thickness/2.,scinti_tile_z_old/2.);

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
  G4double trap_y = params->scinti_tile_thickness*2;
  G4VSolid *cuttrapezoid = new G4Trap("Cuttrapezoid",trap_y,cuttrapezoid_x,cuttrapezoid_z_long+z_offset,cuttrapezoid_z_short+z_offset);
  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateX(-90*deg);
  rotm->rotateZ(90*deg);
  // G4LogicalVolume* checksolid = new G4LogicalVolume(cuttrapezoid,G4Material::GetMaterial("G4_POLYSTYRENE"),"DISPLAYLOGICAL", 0, 0, 0);
  // G4VisAttributes* visattchk = new G4VisAttributes();
  // visattchk->SetVisibility(true);
  // visattchk->SetForceSolid(false);
  // checksolid->SetVisAttributes(visattchk);
  // new G4PVPlacement(rotm,G4ThreeVector(-scinti_tile_x_old/2.,0,scinti_tile_z_old/2.-(cuttrapezoid_z_short+cuttrapezoid_z_long)/4.+z_offset/2.),checksolid,"DISPLAYVOL",hcalenvelope, 0, false, overlapcheck);
  //    visattchk->SetColour(G4Colour::Yellow());

  G4VSolid *scinti_tile_1 = new G4SubtractionSolid("scinti_tile_1",scinti_tile,cuttrapezoid,rotm,G4ThreeVector(-scinti_tile_x_old/2.,0,scinti_tile_z_old/2.-(cuttrapezoid_z_short+cuttrapezoid_z_long)/4.+z_offset/2.));
  rotm = new G4RotationMatrix();
  rotm->rotateX(90*deg);
  rotm->rotateZ(90*deg);
  G4VSolid *scinti_tile_2 = new G4SubtractionSolid("scinti_tile_2",scinti_tile_1,cuttrapezoid,rotm,G4ThreeVector(-scinti_tile_x_old/2.,0,-scinti_tile_z_old/2.+(cuttrapezoid_z_short+cuttrapezoid_z_long)/4.-z_offset/2.));
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

void
PHG4OuterHcalDetector::ConstructHcalSingleScintillators(G4LogicalVolume* hcalenvelope)
{
  G4VSolid *bigtile = ConstructScintillatorBox(hcalenvelope);
  // eta->theta
  G4double delta_eta = params->scinti_eta_coverage / params->n_scinti_tiles;
  G4double eta = 0;
  G4double theta;
  G4double x[4];
  G4double z[4];
  ostringstream name;
  double overhang = (scinti_tile_x - (params->outer_radius - params->inner_radius)) / 2.;
  double offset = 1 * cm + overhang; // add 1cm to make sure the G4ExtrudedSolid
  // is larger than the tile so we do not have
  // funny edge effects when overlapping vols
  double magnet_cutout_x = params->magnet_cutout/cos(params->tilt_angle/rad);
  double x_inner = params->inner_radius - overhang;
  double inner_offset = offset;
  for (int i = 0; i < params->n_scinti_tiles; i++)
    {
      if (i >= params->magnet_cutout_first_scinti)
	{
          x_inner = params->inner_radius - overhang + magnet_cutout_x;
	  inner_offset = offset - magnet_cutout_x; 
	}
      cout << "tile " << i << " starting at " << x_inner/cm << endl;
      theta = M_PI / 2 - PHG4Utils::get_theta(eta); // theta = 90 for eta=0
      x[0] = x_inner;
      z[0] = tan(theta) * params->inner_radius;
      x[1] = params->outer_radius + overhang; // since the tile is tilted, x is not at the outer radius but beyond
      z[1] = tan(theta) * params->outer_radius;
      eta += delta_eta;
      theta = M_PI / 2 - PHG4Utils::get_theta(eta); // theta = 90 for eta=0
      x[2] = x_inner;
      z[2] =  tan(theta) * params->inner_radius;
      x[3] =  params->outer_radius + overhang; // since the tile is tilted, x is not at the outer radius but beyond
      z[3] = tan(theta) * params->outer_radius;
      // apply gap between scintillators
      z[0] += params->scinti_gap_neighbor / 2.;
      z[1] += params->scinti_gap_neighbor / 2.;
      z[2] -= params->scinti_gap_neighbor / 2.;
      z[3] -= params->scinti_gap_neighbor / 2.;
      Point_2 leftsidelow(z[0], x[0]);
      Point_2 leftsidehigh(z[1], x[1]);
      x[0] = params->inner_radius - inner_offset;
      z[0] = x_at_y(leftsidelow, leftsidehigh, x[0]);
      x[1] = params->outer_radius + offset;
      z[1] = x_at_y(leftsidelow, leftsidehigh, x[1]);
      Point_2 rightsidelow(z[2], x[2]);
      Point_2 rightsidehigh(z[3], x[3]);
      x[2] = params->outer_radius + offset;
      z[2] = x_at_y(rightsidelow, rightsidehigh, x[2]);
      x[3] = params->inner_radius - inner_offset;
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
					      params->scinti_tile_thickness + 0.2 * mm,
					      zero, 1.0,
					      zero, 1.0);
      G4RotationMatrix *rotm = new G4RotationMatrix();
      rotm->rotateX(-90 * deg);
      name.str("");
      name << "scintillator_" << i << "_left";
      G4VSolid *scinti_tile =  new G4IntersectionSolid(name.str(), bigtile, scinti, rotm, G4ThreeVector(-(params->inner_radius + params->outer_radius) / 2., 0, 0));
      scinti_tiles_vec[i + params->n_scinti_tiles] = scinti_tile;
      rotm = new G4RotationMatrix();
      rotm->rotateX(90 * deg);
      name.str("");
      name << "scintillator_" << i << "_right";
      scinti_tile =  new G4IntersectionSolid(name.str(), bigtile, scinti, rotm, G4ThreeVector(-(params->inner_radius + params->outer_radius) / 2., 0, 0));
      scinti_tiles_vec[params->n_scinti_tiles - i - 1] =  scinti_tile;
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

G4double
PHG4OuterHcalDetector::x_at_y(Point_2 &p0, Point_2 &p1, G4double yin)
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

int
PHG4OuterHcalDetector::ConstructHcalSingleScintillator_A(G4LogicalVolume* hcalenvelope)
{
  G4VSolid *bigtile = ConstructHcalScintillator(hcalenvelope);
  G4double delta_eta = scinti_eta_coverage/params->n_scinti_tiles;
  G4double offset = 10*cm;
  G4double inner_reference = envelope_inner_radius-offset;
  G4double outer_reference = envelope_outer_radius+offset;
  G4double x[4];
  G4double z[4];
  G4double eta = 0;
  G4double theta;

  for (int j=0; j<params->n_scinti_tiles;j++)
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
      z[0] += params->scinti_gap_neighbor/2.;
      z[1] += params->scinti_gap_neighbor/2.;
      z[2] -= params->scinti_gap_neighbor/2.;
      z[3] -= params->scinti_gap_neighbor/2.;
      vector<G4TwoVector> vertexes;
      for (int i=0; i<4; i++)
	{
	  G4TwoVector v(x[i],z[i]);
	  vertexes.push_back(v);
	}
      G4TwoVector zero(0, 0);
      G4VSolid *scinti =  new G4ExtrudedSolid("ScintillatorTile",
					      vertexes,
					      params->scinti_tile_thickness,
					      zero, 1.0,
					      zero, 1.0);
      G4RotationMatrix *rotm = new G4RotationMatrix();
      rotm->rotateX(-90*deg);
      G4VSolid *scinti_tile =  new G4IntersectionSolid("scintillator",bigtile,scinti,rotm,G4ThreeVector(-(inner_reference+outer_reference)/2., 0, 0));
      scinti_tiles_vec[j+params->n_scinti_tiles] = scinti_tile;
      rotm = new G4RotationMatrix();
      rotm->rotateX(90*deg);
      scinti_tile =  new G4IntersectionSolid("scintillator",bigtile,scinti,rotm,G4ThreeVector(-(inner_reference+outer_reference)/2., 0, 0));
      scinti_tiles_vec[params->n_scinti_tiles-j-1] =  scinti_tile;
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
  ConstructHcalSingleScintillators(hcalenvelope);
  G4AssemblyVolume *assmeblyvol = new G4AssemblyVolume();
  ostringstream name;
  G4ThreeVector g4vec;
  for (unsigned int i=0; i<scinti_tiles_vec.size(); i++)
    {
      name.str("");
      name << scintilogicnameprefix << i;
      G4UserLimits *g4userlimits = NULL;
      if (isfinite(params->steplimits))
	{
	  g4userlimits = new G4UserLimits(params->steplimits);
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
  if (params->active)
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
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeomv3(envelope_inner_radius / cm, (params->place_in_z - steel_plate_z / 2.) / cm, (params->place_in_z + steel_plate_z / 2.) / cm, (envelope_outer_radius-envelope_inner_radius) / cm, n_steel_plates,  params->tilt_angle/rad, 0);
      geo->AddLayerGeom(layer, mygeom);
      if (verbosity > 0) geo->identify();
    }
}

int
PHG4OuterHcalDetector::ConsistencyCheck() const
{
  // just make sure the parameters make a bit of sense
  if (params->inner_radius >= params->outer_radius)
    {
      cout << PHWHERE << ": Inner Radius " << params->inner_radius/cm
	   << " cm larger than Outer Radius " << params->outer_radius/cm
	   << " cm" << endl;
      exit(1);
    }
  if (params->scinti_tile_thickness > params->scinti_gap)
    {
      cout << PHWHERE << "Scintillator thickness " << params->scinti_tile_thickness/cm
	   << " cm larger than scintillator gap " << params->scinti_gap/cm
	   << " cm" << endl;
      exit(1);
    }
  return 0;
}

void
PHG4OuterHcalDetector::SetTiltViaNcross()
{
  if (! params->ncross)
    {
      return;
    }
  if ((isfinite(params->tilt_angle))&&(verbosity > 0))
    {
      cout << "both number of crossings and tilt angle are set" << endl;
      cout << "using number of crossings to determine tilt angle" << endl;
    }
  double mid_radius = params->inner_radius + (params->outer_radius-params->inner_radius)/2.;
  double deltaphi = (2*M_PI/params->n_scinti_plates)*params->ncross;
  Point_2 pnull(0,0);
  Point_2 plow(params->inner_radius,0);
  Point_2 phightmp(1,tan(deltaphi));
  Point_2 pin1(params->inner_radius,0), pin2(0,params->inner_radius),pin3(-params->inner_radius,0);
  Circle_2 inner_circle(pin1,pin2,pin3);
  Point_2 pmid1(mid_radius,0), pmid2(0,mid_radius),pmid3(-mid_radius,0);
  Circle_2 mid_circle(pmid1,pmid2,pmid3);
  Point_2 pout1(params->outer_radius,0), pout2(0,params->outer_radius),pout3(-params->outer_radius,0);
  Circle_2 outer_circle(pout1,pout2,pout3);
  Line_2 l_up(pnull,phightmp);
  vector< CGAL::Object > res;
  CGAL::intersection(outer_circle, l_up, std::back_inserter(res));
  Point_2 upperright;
  vector< CGAL::Object >::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); ++iter)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) >  0)
	    {
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      upperright = pntmp;
	    }
	}
      else
	{
	  cout << "CGAL::Object type not pair..." << endl;
	  exit(1);
	}
    }
  Line_2 l_right(plow,upperright);
  res.clear();
  Point_2 midpoint;
  CGAL::intersection(mid_circle, l_right, std::back_inserter(res));
  for (iter = res.begin(); iter != res.end(); ++iter)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<Circular_k>, unsigned> >(&obj))
	{
	  if (CGAL::to_double(point->first.x()) >  0)
	    {
	      Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
	      midpoint = pntmp;
	    }
	}
      else
	{
	  cout << "CGAL::Object type not pair..." << endl;
	  exit(1);
	}
    }
  // length left side
  double ll = sqrt((CGAL::to_double(midpoint.x()) - params->inner_radius)*(CGAL::to_double(midpoint.x()) - params->inner_radius) + CGAL::to_double(midpoint.y())*CGAL::to_double(midpoint.y()));
  double upside = sqrt(CGAL::to_double(midpoint.x())*CGAL::to_double(midpoint.x()) + CGAL::to_double(midpoint.y())*CGAL::to_double(midpoint.y()));
  //  c^2 = a^2+b^2 - 2ab*cos(gamma)
  // gamma = acos((a^2+b^2=c^2)/2ab
  double tiltangle = acos((ll*ll + upside*upside-params->inner_radius*params->inner_radius)/(2*ll*upside));
  params->tilt_angle = copysign(tiltangle,params->ncross);
  return;
}

// check if tilt angle is reasonable - too large, no intersections with inner radius
int
PHG4OuterHcalDetector::CheckTiltAngle() const
{
  if (fabs(params->tilt_angle) >= M_PI)
    {
      cout << PHWHERE << "invalid tilt angle, abs(tilt) >= 90 deg: " << (params->tilt_angle / deg)
	   << endl;
      exit(1);
    }

  double mid_radius = params->inner_radius + (params->outer_radius - params->inner_radius) / 2.;
  Point_2 pmid(mid_radius, 0); // center of scintillator
  double xcoord = 0;
  double ycoord = mid_radius * tan(params->tilt_angle / rad) ;
  Point_2 pxnull(xcoord, ycoord);
  Line_2 s2(pmid, pxnull);
  Point_2 sc1(params->inner_radius, 0), sc2(0, params->inner_radius), sc3(-params->inner_radius, 0);
  Circle_2 inner_circle(sc1, sc2, sc3);
  vector< CGAL::Object > res;
  CGAL::intersection(inner_circle, s2, std::back_inserter(res));
  if (res.size() == 0)
    {
      cout << PHWHERE << " Tilt angle " << (params->tilt_angle / deg)
	   << " too large, no intersection with inner radius" << endl;
      exit(1);
    }
  return 0;
}
