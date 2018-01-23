#include "PHG4InnerHcalDetector.h"
#include "PHG4HcalDefs.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Circular_kernel_intersections.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Object.h>
#include <CGAL/point_generators_2.h>

#include <cmath>
#include <sstream>

typedef CGAL::Circle_2<PHG4InnerHcalDetector::Circular_k> Circle_2;
typedef CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k> Circular_arc_point_2;
typedef CGAL::Line_2<PHG4InnerHcalDetector::Circular_k> Line_2;
typedef CGAL::Segment_2<PHG4InnerHcalDetector::Circular_k> Segment_2;

using namespace std;

// there is still a minute problem for very low tilt angles where the scintillator
// face touches the boundary instead of the corner, subtracting 1 permille from the total
// scintilator length takes care of this
static double subtract_from_scinti_x = 0.1 * mm;

PHG4InnerHcalDetector::PHG4InnerHcalDetector(PHCompositeNode *Node, PHG4Parameters *parameters, const std::string &dnam)
  : PHG4Detector(Node, dnam)
  , params(parameters)
  , scinti_mother_assembly(nullptr)
  , inner_radius(params->get_double_param("inner_radius") * cm)
  , outer_radius(params->get_double_param("outer_radius") * cm)
  , size_z(params->get_double_param("size_z") * cm)
  , scinti_tile_x(NAN)
  , scinti_tile_x_lower(NAN)
  , scinti_tile_x_upper(NAN)
  , scinti_tile_z(size_z)
  , scinti_tile_thickness(params->get_double_param("scinti_tile_thickness") * cm)
  , scinti_inner_gap(params->get_double_param("scinti_inner_gap") * cm)
  , scinti_outer_gap(params->get_double_param("scinti_outer_gap") * cm)
  , scinti_outer_radius(params->get_double_param("scinti_outer_radius") * cm)
  , tilt_angle(params->get_double_param("tilt_angle") * deg)
  , envelope_inner_radius(inner_radius)
  , envelope_outer_radius(outer_radius)
  , envelope_z(size_z)
  , volume_envelope(NAN)
  , volume_steel(NAN)
  , volume_scintillator(NAN)
  , n_scinti_plates(params->get_int_param(PHG4HcalDefs::scipertwr) * params->get_int_param("n_towers"))
  , n_scinti_tiles(params->get_int_param("n_scinti_tiles"))
  , active(params->get_int_param("active"))
  , absorberactive(params->get_int_param("absorberactive"))
  , layer(0)
  , scintilogicnameprefix("HcalInnerScinti")
{
  // allocate memory for scintillator plates
  scinti_tiles_vec.assign(2 * n_scinti_tiles, static_cast<G4VSolid *>(nullptr));
}

PHG4InnerHcalDetector::~PHG4InnerHcalDetector()
{
  delete scinti_mother_assembly;
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4InnerHcalDetector::IsInInnerHcal(G4VPhysicalVolume *volume) const
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

G4VSolid *
PHG4InnerHcalDetector::ConstructScintillatorBox(G4LogicalVolume *hcalenvelope)
{
  double mid_radius = inner_radius + (outer_radius - inner_radius) / 2.;
  Point_2 p_in_1(mid_radius, 0);  // center of scintillator

  // length of upper edge (middle till outer circle intersect
  // x/y coordinate of end of center vertical
  double xcoord = scinti_tile_thickness / 2. * sin(fabs(tilt_angle) / rad) + mid_radius;
  double ycoord = scinti_tile_thickness / 2. * cos(fabs(tilt_angle) / rad) + 0;
  Point_2 p_upperedge(xcoord, ycoord);
  Line_2 s2(p_in_1, p_upperedge);  // center vertical

  Line_2 perp = s2.perpendicular(p_upperedge);
  Point_2 sc1(outer_radius, 0), sc2(0, outer_radius), sc3(-outer_radius, 0);
  Circle_2 outer_circle(sc1, sc2, sc3);
  vector<CGAL::Object> res;
  CGAL::intersection(outer_circle, perp, std::back_inserter(res));
  Point_2 upperright;
  vector<CGAL::Object>::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> >(&obj))
    {
      if (CGAL::to_double(point->first.x()) > CGAL::to_double(p_upperedge.x()))
      {
        double deltax = CGAL::to_double(point->first.x()) - CGAL::to_double(p_upperedge.x());
        double deltay = CGAL::to_double(point->first.y()) - CGAL::to_double(p_upperedge.y());
        // the scintillator is twice as long
        scinti_tile_x_upper = sqrt(deltax * deltax + deltay * deltay);  //
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
  xcoord = mid_radius - scinti_tile_thickness / 2. * sin(fabs(tilt_angle) / rad);
  ycoord = 0 - scinti_tile_thickness / 2. * cos(fabs(tilt_angle) / rad);
  Point_2 p_loweredge(xcoord, ycoord);
  Line_2 s3(p_in_1, p_loweredge);
  Line_2 l_lower = s3.perpendicular(p_loweredge);
  Point_2 ic1(inner_radius, 0), ic2(0, inner_radius), ic3(-inner_radius, 0);
  Circle_2 inner_circle(ic1, ic2, ic3);
  res.clear();
  CGAL::intersection(inner_circle, l_lower, std::back_inserter(res));
  Point_2 lowerleft;
  // we have 2 intersections - we want the one furthest to the right (largest x). The correct one is
  // certainly > 0 but the second one depends on the tilt angle and might also be > 0
  double minx = 0;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> >(&obj))
    {
      if (CGAL::to_double(point->first.x()) > minx)
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
  scinti_tile_x = scinti_tile_x_upper + scinti_tile_x_lower - ((outer_radius - scinti_outer_radius) / cos(tilt_angle / rad));
  scinti_tile_x -= subtract_from_scinti_x;
  G4VSolid *scintibox = new G4Box("ScintiTile", scinti_tile_x / 2., scinti_tile_thickness / 2., scinti_tile_z / 2.);
  volume_scintillator = scintibox->GetCubicVolume() * n_scinti_plates;
  return scintibox;
}

G4VSolid *
PHG4InnerHcalDetector::ConstructSteelPlate(G4LogicalVolume *hcalenvelope)
{
  // calculate steel plate on top of the scinti box. Lower edge is the upper edge of
  // the scintibox + 1/2 the airgap
  double mid_radius = inner_radius + (outer_radius - inner_radius) / 2.;
  // first the lower edge, just like the scinti box, just add the air gap
  // and calculate intersection of edge with inner and outer radius.
  Point_2 p_in_1(mid_radius, 0);  // center of lower scintillator
  double angle_mid_scinti = M_PI / 2. + tilt_angle / rad;
  double xcoord = scinti_inner_gap / 2. * cos(angle_mid_scinti / rad) + mid_radius;
  double ycoord = scinti_inner_gap / 2. * sin(angle_mid_scinti / rad) + 0;
  Point_2 p_loweredge(xcoord, ycoord);
  Line_2 s2(p_in_1, p_loweredge);               // center vertical
  Line_2 perp = s2.perpendicular(p_loweredge);  // that is the lower edge of the steel plate
  Point_2 sc1(inner_radius, 0), sc2(0, inner_radius), sc3(-inner_radius, 0);
  Circle_2 inner_circle(sc1, sc2, sc3);
  vector<CGAL::Object> res;
  CGAL::intersection(inner_circle, perp, std::back_inserter(res));
  Point_2 lowerleft;
  vector<CGAL::Object>::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> >(&obj))
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

  double xcoord2 = scinti_outer_gap / 2. * cos(angle_mid_scinti / rad) + mid_radius;
  double ycoord2 = scinti_outer_gap / 2. * sin(angle_mid_scinti / rad) + 0;
  Point_2 p_loweredge2(xcoord2, ycoord2);
  Line_2 s2_2(p_in_1, p_loweredge2);                // center vertical
  Line_2 perp2 = s2_2.perpendicular(p_loweredge2);  // that is the lower edge of the steel plate

  Point_2 so1(outer_radius, 0), so2(0, outer_radius), so3(-outer_radius, 0);
  Circle_2 outer_circle(so1, so2, so3);
  res.clear();  // just clear the content from the last intersection search
  CGAL::intersection(outer_circle, perp2, std::back_inserter(res));
  Point_2 lowerright;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> >(&obj))
    {
      if (CGAL::to_double(point->first.x()) > CGAL::to_double(p_loweredge.x()))
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

  double phi_midpoint = 2 * M_PI / n_scinti_plates;
  double xmidpoint = cos(phi_midpoint) * mid_radius;
  double ymidpoint = sin(phi_midpoint) * mid_radius;
  // angle of perp line at center of scintillator
  angle_mid_scinti = (M_PI / 2. - phi_midpoint) - (M_PI / 2. + tilt_angle / rad);
  double xcoordup = xmidpoint - scinti_inner_gap / 2. * sin(angle_mid_scinti / rad);
  double ycoordup = ymidpoint - scinti_inner_gap / 2. * cos(angle_mid_scinti / rad);
  Point_2 upperleft;
  Point_2 upperright;
  Point_2 mid_upperscint(xmidpoint, ymidpoint);
  Point_2 p_upperedge(xcoordup, ycoordup);
  {
    Line_2 sup(mid_upperscint, p_upperedge);       // center vertical
    Line_2 perp = sup.perpendicular(p_upperedge);  // that is the upper edge of the steel plate
    Point_2 sc1(inner_radius, 0), sc2(0, inner_radius), sc3(-inner_radius, 0);
    Circle_2 inner_circle(sc1, sc2, sc3);
    vector<CGAL::Object> res;
    CGAL::intersection(inner_circle, perp, std::back_inserter(res));
    vector<CGAL::Object>::const_iterator iter;
    double pxmax = 0.;
    for (iter = res.begin(); iter != res.end(); ++iter)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> >(&obj))
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

    double xcoordup2 = xmidpoint - scinti_outer_gap / 2. * sin(angle_mid_scinti / rad);
    double ycoordup2 = ymidpoint - scinti_outer_gap / 2. * cos(angle_mid_scinti / rad);
    Point_2 p_upperedge2(xcoordup2, ycoordup2);
    Line_2 sup2(mid_upperscint, p_upperedge2);        // center vertical
    Line_2 perp2 = sup2.perpendicular(p_upperedge2);  // that is the upper edge of the steel plate

    Point_2 so1(outer_radius, 0), so2(0, outer_radius), so3(-outer_radius, 0);
    Circle_2 outer_circle(so1, so2, so3);
    res.clear();  // just clear the content from the last intersection search
    CGAL::intersection(outer_circle, perp2, std::back_inserter(res));
    for (iter = res.begin(); iter != res.end(); ++iter)
    {
      CGAL::Object obj = *iter;
      if (const std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> >(&obj))
      {
        if (CGAL::to_double(point->first.x()) > CGAL::to_double(p_loweredge.x()))
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
  G4VSolid *steel_plate = new G4ExtrudedSolid("SteelPlate",
                                              vertexes,
                                              size_z / 2.0,
                                              zero, 1.0,
                                              zero, 1.0);

  //  DisplayVolume(steel_plate, hcalenvelope);
  volume_steel = steel_plate->GetCubicVolume() * n_scinti_plates;
  return steel_plate;
}

void PHG4InnerHcalDetector::ShiftSecantToTangent(Point_2 &lowleft, Point_2 &upleft, Point_2 &upright, Point_2 &lowright)
{
  Line_2 secant(lowleft, upleft);
  Segment_2 upedge(upleft, upright);
  Segment_2 lowedge(lowleft, lowright);
  double xmid = (CGAL::to_double(lowleft.x()) + CGAL::to_double(upleft.x())) / 2.;
  double ymid = (CGAL::to_double(lowleft.y()) + CGAL::to_double(upleft.y())) / 2.;
  Point_2 midpoint(xmid, ymid);
  Line_2 sekperp = secant.perpendicular(midpoint);
  Point_2 sc1(inner_radius, 0), sc2(0, inner_radius), sc3(-inner_radius, 0);
  Circle_2 inner_circle(sc1, sc2, sc3);
  vector<CGAL::Object> res;
  CGAL::intersection(inner_circle, sekperp, std::back_inserter(res));
  vector<CGAL::Object>::const_iterator iter;
  double pxmax = 0.;
  Point_2 tangtouch;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> >(&obj))
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

// Construct the envelope and the call the
// actual inner hcal construction
void PHG4InnerHcalDetector::Construct(G4LogicalVolume *logicWorld)
{
  G4Material *Air = G4Material::GetMaterial("G4_AIR");
  G4VSolid *hcal_envelope_cylinder = new G4Tubs("InnerHcal_envelope_solid", envelope_inner_radius, envelope_outer_radius, envelope_z / 2., 0, 2 * M_PI);
  volume_envelope = hcal_envelope_cylinder->GetCubicVolume();
  G4LogicalVolume *hcal_envelope_log = new G4LogicalVolume(hcal_envelope_cylinder, Air, G4String("Hcal_envelope"), 0, 0, 0);
  G4VisAttributes *hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(false);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::White());
  hcal_envelope_log->SetVisAttributes(hcalVisAtt);
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(params->get_double_param("rot_x") * deg);
  hcal_rotm.rotateY(params->get_double_param("rot_y") * deg);
  hcal_rotm.rotateZ(params->get_double_param("rot_z") * deg);
  new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(params->get_double_param("place_x") * cm, params->get_double_param("place_y") * cm, params->get_double_param("place_z") * cm)), hcal_envelope_log, "InnerHcalEnvelope", logicWorld, 0, false, overlapcheck);
  ConstructInnerHcal(hcal_envelope_log);
  return;
}

int PHG4InnerHcalDetector::ConstructInnerHcal(G4LogicalVolume *hcalenvelope)
{
  ConsistencyCheck();
  SetTiltViaNcross();  // if number of crossings is set, use it to determine tilt
  CheckTiltAngle();    // die if the tilt angle is out of range
  G4VSolid *steel_plate = ConstructSteelPlate(hcalenvelope);
  G4LogicalVolume *steel_logical = new G4LogicalVolume(steel_plate, G4Material::GetMaterial(params->get_string_param("material")), "HcalInnerSteelPlate", 0, 0, 0);
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(true);
  visattchk->SetColour(G4Colour::Grey());
  steel_logical->SetVisAttributes(visattchk);
  scinti_mother_assembly = ConstructHcalScintillatorAssembly(hcalenvelope);
  double phi = 0;
  double deltaphi = 2 * M_PI / n_scinti_plates;
  ostringstream name;
  double middlerad = outer_radius - (outer_radius - inner_radius) / 2.;
  // for shorter scintillators we have to shift by 1/2 the "missing length"
  double scinti_tile_orig_length = scinti_tile_x_upper + scinti_tile_x_lower - subtract_from_scinti_x;
  double shiftslat = fabs(scinti_tile_x_lower - scinti_tile_x_upper) / 2. + (scinti_tile_orig_length - scinti_tile_x) / 2.;
  // calculate phi offset (copied from code inside following loop):
  // first get the center point (phi=0) so it's middlerad/0
  // then shift the scintillator center as documented in loop
  // then
  // for positive tilt angles we need the lower left corner of the scintillator
  // for negative tilt angles we need the upper right corner of the scintillator
  // as it turns out the code uses the middle of the face of the scintillator
  // as reference, if this is a problem the code needs to be modified to
  // actually calculate the corner (but the math of the construction is that
  // the middle of the scintillator sits at zero)
  double xp = cos(phi) * middlerad;
  double yp = sin(phi) * middlerad;
  xp -= cos((-tilt_angle) / rad - phi) * shiftslat;
  yp += sin((-tilt_angle) / rad - phi) * shiftslat;
  if (tilt_angle > 0)
  {
    double xo = xp - (scinti_tile_x / 2.) * cos(tilt_angle / rad);
    double yo = yp - (scinti_tile_x / 2.) * sin(tilt_angle / rad);
    phi = -atan(yo / xo);
  }
  else if (tilt_angle < 0)
  {
    double xo = xp + (scinti_tile_x / 2.) * cos(tilt_angle / rad);
    double yo = yp + (scinti_tile_x / 2.) * sin(tilt_angle / rad);
    phi = -atan(yo / xo);
  }
  // else (for tilt_angle = 0) phi stays zero
  for (int i = 0; i < n_scinti_plates; i++)
  {
    G4RotationMatrix *Rot = new G4RotationMatrix();
    double ypos = sin(phi) * middlerad;
    double xpos = cos(phi) * middlerad;
    // the center of the scintillator is not the center of the inner hcal
    // but depends on the tilt angle. Therefore we need to shift
    // the center from the mid point
    ypos += sin((-tilt_angle) / rad - phi) * shiftslat;
    xpos -= cos((-tilt_angle) / rad - phi) * shiftslat;
    Rot->rotateZ(phi * rad + tilt_angle);
    G4ThreeVector g4vec(xpos, ypos, 0);
    // great the MakeImprint always adds 1 to the copy number and 0 has a
    // special meaning (which then also adds 1). Basically our volume names
    // will start at 1 instead of 0 and there is nothing short of patching this
    // method. I'll take care of this in the decoding of the volume name
    // AAAAAAARHGS
    scinti_mother_assembly->MakeImprint(hcalenvelope, g4vec, Rot, i, overlapcheck);
    Rot = new G4RotationMatrix();
    Rot->rotateZ(-phi * rad);
    name.str("");
    name << "InnerHcalSteel_" << i;
    steel_absorber_vec.insert(new G4PVPlacement(Rot, G4ThreeVector(0, 0, 0), steel_logical, name.str().c_str(), hcalenvelope, 0, i, overlapcheck));
    phi += deltaphi;
  }
  return 0;
}

// split the big scintillator into tiles covering eta ranges
// since they are tilted it is not a straightforward theta cut
// it starts at a given eta at the inner radius but the outer radius needs adjusting
void PHG4InnerHcalDetector::ConstructHcalSingleScintillators(G4LogicalVolume *hcalenvelope)
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
  double offset = 1 * cm + overhang;  // add 1cm to make sure the G4ExtrudedSolid
  // is larger than the tile so we do not have
  // funny edge effects when overlapping vols
  double scinti_gap_neighbor = params->get_double_param("scinti_gap_neighbor") * cm;
  for (int i = 0; i < n_scinti_tiles; i++)
  {
    theta = M_PI / 2 - PHG4Utils::get_theta(eta);  // theta = 90 for eta=0
    x[0] = inner_radius - overhang;
    z[0] = tan(theta) * inner_radius;
    x[1] = outer_radius + overhang;  // since the tile is tilted, x is not at the outer radius but beyond
    z[1] = tan(theta) * outer_radius;
    eta += delta_eta;
    theta = M_PI / 2 - PHG4Utils::get_theta(eta);  // theta = 90 for eta=0
    x[2] = inner_radius - overhang;
    z[2] = tan(theta) * inner_radius;
    x[3] = outer_radius + overhang;  // since the tile is tilted, x is not at the outer radius but beyond
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

    G4VSolid *scinti = new G4ExtrudedSolid("ScintillatorTile",
                                           vertexes,
                                           scinti_tile_thickness + 0.2 * mm,
                                           zero, 1.0,
                                           zero, 1.0);
    G4RotationMatrix *rotm = new G4RotationMatrix();
    rotm->rotateX(-90 * deg);
    name.str("");
    name << "scintillator_" << i << "_left";
    G4VSolid *scinti_tile = new G4IntersectionSolid(name.str(), bigtile, scinti, rotm, G4ThreeVector(-(inner_radius + outer_radius) / 2., 0, 0));
    scinti_tiles_vec[i + n_scinti_tiles] = scinti_tile;
    rotm = new G4RotationMatrix();
    rotm->rotateX(90 * deg);
    name.str("");
    name << "scintillator_" << i << "_right";
    scinti_tile = new G4IntersectionSolid(name.str(), bigtile, scinti, rotm, G4ThreeVector(-(inner_radius + outer_radius) / 2., 0, 0));
    scinti_tiles_vec[n_scinti_tiles - i - 1] = scinti_tile;
  }

  // for (unsigned int i=0; i<scinti_tiles_vec.size(); i++)
  //    {
  //      if (scinti_tiles_vec[i])
  //  	 {
  //   	   DisplayVolume(scinti_tiles_vec[i],hcalenvelope );
  //   	 }
  //     }

  return;
}

double
PHG4InnerHcalDetector::x_at_y(Point_2 &p0, Point_2 &p1, double yin)
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
  if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
  {
    xret = CGAL::to_double(ipoint->x());
  }
  else
  {
    cout << PHWHERE << " failed for y = " << y << endl;
    cout << "p0(x): " << CGAL::to_double(p0.x()) << ", p0(y): " << CGAL::to_double(p0.y()) << endl;
    cout << "p1(x): " << CGAL::to_double(p1.x()) << ", p1(y): " << CGAL::to_double(p1.y()) << endl;
    exit(1);
  }
  return xret;
}

G4AssemblyVolume *
PHG4InnerHcalDetector::ConstructHcalScintillatorAssembly(G4LogicalVolume *hcalenvelope)
{
  ConstructHcalSingleScintillators(hcalenvelope);
  G4AssemblyVolume *assmeblyvol = new G4AssemblyVolume();
  ostringstream name;
  G4ThreeVector g4vec;

  double steplimits = params->get_double_param("steplimits") * cm;
  for (unsigned int i = 0; i < scinti_tiles_vec.size(); i++)
  {
    name.str("");
    name << scintilogicnameprefix << i;
    G4UserLimits *g4userlimits = nullptr;
    if (isfinite(steplimits))
    {
      g4userlimits = new G4UserLimits(steplimits);
    }
    G4LogicalVolume *scinti_tile_logic = new G4LogicalVolume(scinti_tiles_vec[i], G4Material::GetMaterial("G4_POLYSTYRENE"), name.str().c_str(), nullptr, nullptr, g4userlimits);
    G4VisAttributes *visattchk = new G4VisAttributes();
    visattchk->SetVisibility(true);
    visattchk->SetForceSolid(true);
    visattchk->SetColour(G4Colour::Green());
    scinti_tile_logic->SetVisAttributes(visattchk);
    assmeblyvol->AddPlacedVolume(scinti_tile_logic, g4vec, nullptr);
  }
  return assmeblyvol;
}

int PHG4InnerHcalDetector::DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm)
{
  static int i = 0;
  G4LogicalVolume *checksolid = new G4LogicalVolume(volume, G4Material::GetMaterial("G4_POLYSTYRENE"), "DISPLAYLOGICAL", 0, 0, 0);
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  switch (i)
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
// check if tilt angle is reasonable - too large, no intersections with inner radius
int PHG4InnerHcalDetector::CheckTiltAngle() const
{
  if (fabs(tilt_angle / rad) >= M_PI)
  {
    cout << PHWHERE << "invalid tilt angle, abs(tilt) >= 90 deg: " << (tilt_angle / deg)
         << endl;
    exit(1);
  }

  double mid_radius = inner_radius + (outer_radius - inner_radius) / 2.;
  Point_2 pmid(mid_radius, 0);  // center of scintillator
  double xcoord = 0;
  double ycoord = mid_radius * tan(tilt_angle / rad);
  Point_2 pxnull(xcoord, ycoord);
  Line_2 s2(pmid, pxnull);
  Point_2 sc1(inner_radius, 0), sc2(0, inner_radius), sc3(-inner_radius, 0);
  Circle_2 inner_circle(sc1, sc2, sc3);
  vector<CGAL::Object> res;
  CGAL::intersection(inner_circle, s2, std::back_inserter(res));
  if (res.size() == 0)
  {
    cout << PHWHERE << " Tilt angle " << (tilt_angle / deg)
         << " too large, no intersection with inner radius" << endl;
    exit(1);
  }
  return 0;
}

int PHG4InnerHcalDetector::ConsistencyCheck() const
{
  // just make sure the parameters make a bit of sense
  if (inner_radius >= outer_radius)
  {
    cout << PHWHERE << ": Inner Radius " << inner_radius / cm
         << " cm larger than Outer Radius " << outer_radius / cm
         << " cm" << endl;
    gSystem->Exit(1);
  }
  if (scinti_tile_thickness > scinti_inner_gap)
  {
    cout << PHWHERE << "Scintillator thickness " << scinti_tile_thickness / cm
         << " cm larger than scintillator inner gap " << scinti_inner_gap / cm
         << " cm" << endl;
    gSystem->Exit(1);
  }
  if (scinti_outer_radius <= inner_radius)
  {
    cout << PHWHERE << "Scintillator outer radius " << scinti_outer_radius / cm
         << " cm smaller than inner radius " << inner_radius / cm
         << " cm" << endl;
    gSystem->Exit(1);
  }
  return 0;
}

void PHG4InnerHcalDetector::SetTiltViaNcross()
{
  // if tilt angle is set it takes precedence
  // this happens during readback where the tilt angle was set
  // in this method
  int ncross = params->get_int_param("ncross");
  if (!ncross || isfinite(tilt_angle))
  {
    // flag ncrossing as not used
    params->set_int_param("ncross",0);
    return;
  }
  if ((isfinite(tilt_angle)) && (verbosity > 0))
  {
    cout << "both number of crossings and tilt angle are set" << endl;
    cout << "using number of crossings to determine tilt angle" << endl;
  }
  double mid_radius = inner_radius + (outer_radius - inner_radius) / 2.;
  double deltaphi = (2 * M_PI / n_scinti_plates) * ncross;
  Point_2 pnull(0, 0);
  Point_2 plow(inner_radius, 0);
  Point_2 phightmp(1, tan(deltaphi));
  Point_2 pin1(inner_radius, 0), pin2(0, inner_radius), pin3(-inner_radius, 0);
  Circle_2 inner_circle(pin1, pin2, pin3);
  Point_2 pmid1(mid_radius, 0), pmid2(0, mid_radius), pmid3(-mid_radius, 0);
  Circle_2 mid_circle(pmid1, pmid2, pmid3);
  Point_2 pout1(outer_radius, 0), pout2(0, outer_radius), pout3(-outer_radius, 0);
  Circle_2 outer_circle(pout1, pout2, pout3);
  Line_2 l_up(pnull, phightmp);
  vector<CGAL::Object> res;
  CGAL::intersection(outer_circle, l_up, std::back_inserter(res));
  Point_2 upperright;
  vector<CGAL::Object>::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> >(&obj))
    {
      if (CGAL::to_double(point->first.x()) > 0)
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
  Line_2 l_right(plow, upperright);
  res.clear();
  Point_2 midpoint;
  CGAL::intersection(mid_circle, l_right, std::back_inserter(res));
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4InnerHcalDetector::Circular_k>, unsigned> >(&obj))
    {
      if (CGAL::to_double(point->first.x()) > 0)
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
  double ll = sqrt((CGAL::to_double(midpoint.x()) - inner_radius) * (CGAL::to_double(midpoint.x()) - inner_radius) + CGAL::to_double(midpoint.y()) * CGAL::to_double(midpoint.y()));
  double upside = sqrt(CGAL::to_double(midpoint.x()) * CGAL::to_double(midpoint.x()) + CGAL::to_double(midpoint.y()) * CGAL::to_double(midpoint.y()));
  //  c^2 = a^2+b^2 - 2ab*cos(gamma)
  // gamma = acos((a^2+b^2=c^2)/2ab
  double tiltangle = acos((ll * ll + upside * upside - inner_radius * inner_radius) / (2 * ll * upside));
  tiltangle = tiltangle * rad;
  tilt_angle = copysign(tiltangle, ncross);
  params->set_double_param("tilt_angle", tilt_angle / deg);
  return;
}

void PHG4InnerHcalDetector::Print(const string &what) const
{
  cout << "Inner Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Volume Envelope: " << volume_envelope / cm / cm / cm << " cm^3" << endl;
    cout << "Volume Steel: " << volume_steel / cm / cm / cm << " cm^3" << endl;
    cout << "Volume Scintillator: " << volume_scintillator / cm / cm / cm << " cm^3" << endl;
    cout << "Volume Air: " << (volume_envelope - volume_steel - volume_scintillator) / cm / cm / cm << " cm^3" << endl;
  }
  return;
}
