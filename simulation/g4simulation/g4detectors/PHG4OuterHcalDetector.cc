#include "PHG4OuterHcalDetector.h"

#include "PHG4HcalDefs.h"
#include "PHG4OuterHcalDisplayAction.h"
#include "PHG4OuterHcalFieldSetup.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Subsystem.h>
#include <g4main/PHG4Utils.h>

#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TSystem.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VSolid.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpedantic"
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Circular_kernel_intersections.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Object.h>
#include <CGAL/point_generators_2.h>
#pragma GCC diagnostic pop

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <sstream>

class PHCompositeNode;

using Circle_2 = CGAL::Circle_2<PHG4OuterHcalDetector::Circular_k>;
using Circular_arc_point_2 = CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>;
using Line_2 = CGAL::Line_2<PHG4OuterHcalDetector::Circular_k>;
using Segment_2 = CGAL::Segment_2<PHG4OuterHcalDetector::Circular_k>;

// just for debugging if you want a single layer of scintillators at the center of the world
//#define SCINTITEST

// face touches the boundary instead of the corner, subtracting 1 permille from the total
// scintilator length takes care of this
static double subtract_from_scinti_x = 0.1 * mm;

PHG4OuterHcalDetector::PHG4OuterHcalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parames, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4OuterHcalDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(parames)
  , m_InnerRadius(m_Params->get_double_param("inner_radius") * cm)
  , m_OuterRadius(m_Params->get_double_param("outer_radius") * cm)
  , m_SizeZ(m_Params->get_double_param("size_z") * cm)
  , m_ScintiTileZ(m_SizeZ)
  , m_ScintiTileThickness(m_Params->get_double_param("scinti_tile_thickness") * cm)
  , m_ScintiGap(m_Params->get_double_param("scinti_gap") * cm)
  , m_ScintiInnerRadius(m_Params->get_double_param("scinti_inner_radius") * cm)
  , m_ScintiOuterRadius(m_Params->get_double_param("scinti_outer_radius") * cm)
  , m_TiltAngle(m_Params->get_double_param("tilt_angle") * deg)
  , m_EnvelopeInnerRadius(m_InnerRadius)
  , m_EnvelopeOuterRadius(m_OuterRadius)
  , m_EnvelopeZ(m_SizeZ)
  , m_NumScintiPlates(m_Params->get_int_param(PHG4HcalDefs::scipertwr) * m_Params->get_int_param("n_towers"))
  , m_NumScintiTiles(m_Params->get_int_param("n_scinti_tiles"))
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_ScintiLogicNamePrefix("HcalOuterScinti")
{
  m_ScintiTilesVec.assign(2 * m_NumScintiTiles, static_cast<G4VSolid *>(nullptr));
}

PHG4OuterHcalDetector::~PHG4OuterHcalDetector()
{
  delete m_ScintiMotherAssembly;
  delete m_FieldSetup;
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4OuterHcalDetector::IsInOuterHcal(G4VPhysicalVolume *volume) const
{
  if (m_AbsorberActiveFlag)
  {
    if (m_SteelAbsorberVec.find(volume) != m_SteelAbsorberVec.end())
    {
      return -1;
    }
  }
  if (m_ActiveFlag)
  {
    if (m_ScintiTilePhysVolMap.find(volume) != m_ScintiTilePhysVolMap.end())
    {
      return 1;
    }
  }
  return 0;
}

G4VSolid *
PHG4OuterHcalDetector::ConstructScintillatorBox(G4LogicalVolume * /*hcalenvelope*/)
{
  double mid_radius = m_InnerRadius + (m_OuterRadius - m_InnerRadius) / 2.;
  PHG4OuterHcalDetector::Point_2 p_in_1(mid_radius, 0);  // center of scintillator

  // length of upper edge (middle till outer circle intersect
  // x/y coordinate of end of center vertical
  double xcoord = m_ScintiTileThickness / 2. * sin(fabs(m_TiltAngle) / rad) + mid_radius;
  double ycoord = m_ScintiTileThickness / 2. * cos(fabs(m_TiltAngle) / rad) + 0;
  PHG4OuterHcalDetector::Point_2 p_upperedge(xcoord, ycoord);
  Line_2 s2(p_in_1, p_upperedge);  // center vertical

  Line_2 perp = s2.perpendicular(p_upperedge);
  PHG4OuterHcalDetector::Point_2 sc1(m_OuterRadius, 0), sc2(0, m_OuterRadius), sc3(-m_OuterRadius, 0);
  Circle_2 outer_circle(sc1, sc2, sc3);
  std::vector<CGAL::Object> res;
  CGAL::intersection(outer_circle, perp, std::back_inserter(res));
  PHG4OuterHcalDetector::Point_2 upperright;
  std::vector<CGAL::Object>::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned>>(&obj))
    {
      if (CGAL::to_double(point->first.x()) > CGAL::to_double(p_upperedge.x()))
      {
        double deltax = CGAL::to_double(point->first.x()) - CGAL::to_double(p_upperedge.x());
        double deltay = CGAL::to_double(point->first.y()) - CGAL::to_double(p_upperedge.y());
        // the scintillator is twice as long
        m_ScintiTileXUpper = sqrt(deltax * deltax + deltay * deltay);  //
        PHG4OuterHcalDetector::Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
        upperright = pntmp;
      }
    }
    else
    {
      std::cout << "CGAL::Object type not pair..." << std::endl;
    }
  }
  // length of lower edge (middle till inner circle intersect
  xcoord = mid_radius - m_ScintiTileThickness / 2. * sin(fabs(m_TiltAngle) / rad);
  ycoord = 0 - m_ScintiTileThickness / 2. * cos(fabs(m_TiltAngle) / rad);
  PHG4OuterHcalDetector::Point_2 p_loweredge(xcoord, ycoord);
  Line_2 s3(p_in_1, p_loweredge);
  Line_2 l_lower = s3.perpendicular(p_loweredge);
  PHG4OuterHcalDetector::Point_2 ic1(m_InnerRadius, 0), ic2(0, m_InnerRadius), ic3(-m_InnerRadius, 0);
  Circle_2 inner_circle(ic1, ic2, ic3);
  res.clear();
  CGAL::intersection(inner_circle, l_lower, std::back_inserter(res));
  PHG4OuterHcalDetector::Point_2 lowerleft;
  // we have 2 intersections - we want the one furthest to the right (largest x). The correct one is
  // certainly > 0 but the second one depends on the tilt angle and might also be > 0
  double minx = 0;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned>>(&obj))
    {
      if (CGAL::to_double(point->first.x()) > minx)
      {
        minx = CGAL::to_double(point->first.x());
        double deltax = CGAL::to_double(point->first.x()) - CGAL::to_double(p_loweredge.x());
        double deltay = CGAL::to_double(point->first.y()) - CGAL::to_double(p_loweredge.y());
        m_ScintiTileXLower = sqrt(deltax * deltax + deltay * deltay);
        PHG4OuterHcalDetector::Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
        lowerleft = pntmp;
      }
    }
  }
  m_ScintiTileX = m_ScintiTileXUpper + m_ScintiTileXLower - (m_OuterRadius - m_ScintiOuterRadius) / cos(m_TiltAngle / rad) - (m_ScintiInnerRadius - m_InnerRadius) / cos(m_TiltAngle / rad);
  m_ScintiTileX -= subtract_from_scinti_x;
  G4VSolid *scintibox = new G4Box("ScintiTile", m_ScintiTileX / 2., m_ScintiTileThickness / 2., m_ScintiTileZ / 2.);
  m_VolumeScintillator = scintibox->GetCubicVolume() * m_NumScintiPlates;
  return scintibox;
}

G4VSolid *
PHG4OuterHcalDetector::ConstructSteelPlate(G4LogicalVolume * /*hcalenvelope*/)
{
  // calculate steel plate on top of the scinti box. Lower edge is the upper edge of
  // the scintibox + 1/2 the airgap
  double mid_radius = m_InnerRadius + (m_OuterRadius - m_InnerRadius) / 2.;
  // first the lower edge, just like the scinti box, just add the air gap
  // and calculate intersection of edge with inner and outer radius.
  PHG4OuterHcalDetector::Point_2 p_in_1(mid_radius, 0);  // center of lower scintillator
  double angle_mid_scinti = M_PI / 2. + m_TiltAngle / rad;
  double xcoord = m_ScintiGap / 2. * cos(angle_mid_scinti / rad) + mid_radius;
  double ycoord = m_ScintiGap / 2. * sin(angle_mid_scinti / rad) + 0;
  PHG4OuterHcalDetector::Point_2 p_loweredge(xcoord, ycoord);
  Line_2 s2(p_in_1, p_loweredge);               // center vertical
  Line_2 perp = s2.perpendicular(p_loweredge);  // that is the lower edge of the steel plate
  PHG4OuterHcalDetector::Point_2 sc1(m_InnerRadius, 0), sc2(0, m_InnerRadius), sc3(-m_InnerRadius, 0);
  Circle_2 inner_circle(sc1, sc2, sc3);
  std::vector<CGAL::Object> res;
  CGAL::intersection(inner_circle, perp, std::back_inserter(res));
  PHG4OuterHcalDetector::Point_2 lowerleft;
  std::vector<CGAL::Object>::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned>>(&obj))
    {
      if (CGAL::to_double(point->first.x()) > 0)
      {
        PHG4OuterHcalDetector::Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
        lowerleft = pntmp;
      }
    }
    else
    {
      std::cout << "CGAL::Object type not pair..." << std::endl;
    }
  }
  PHG4OuterHcalDetector::Point_2 so1(m_OuterRadius, 0), so2(0, m_OuterRadius), so3(-m_OuterRadius, 0);
  Circle_2 outer_circle(so1, so2, so3);
  res.clear();  // just clear the content from the last intersection search
  CGAL::intersection(outer_circle, perp, std::back_inserter(res));
  PHG4OuterHcalDetector::Point_2 lowerright;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned>>(&obj))
    {
      if (CGAL::to_double(point->first.x()) > CGAL::to_double(p_loweredge.x()))
      {
        PHG4OuterHcalDetector::Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
        lowerright = pntmp;
      }
    }
    else
    {
      std::cout << "CGAL::Object type not pair..." << std::endl;
    }
  }
  // now we have the lower left and rigth corner, now find the upper edge
  // find the center of the upper scintilator

  double phi_midpoint = 2 * M_PI / m_NumScintiPlates;
  double xmidpoint = cos(phi_midpoint) * mid_radius;
  double ymidpoint = sin(phi_midpoint) * mid_radius;
  // angle of perp line at center of scintillator
  angle_mid_scinti = (M_PI / 2. - phi_midpoint) - (M_PI / 2. + m_TiltAngle / rad);
  double xcoordup = xmidpoint - m_ScintiGap / 2. * sin(angle_mid_scinti / rad);
  double ycoordup = ymidpoint - m_ScintiGap / 2. * cos(angle_mid_scinti / rad);
  PHG4OuterHcalDetector::Point_2 upperleft;
  PHG4OuterHcalDetector::Point_2 upperright;
  PHG4OuterHcalDetector::Point_2 mid_upperscint(xmidpoint, ymidpoint);
  PHG4OuterHcalDetector::Point_2 p_upperedge(xcoordup, ycoordup);
  Line_2 sup(mid_upperscint, p_upperedge);        // center vertical
  Line_2 perpA = sup.perpendicular(p_upperedge);  // that is the upper edge of the steel plate
  PHG4OuterHcalDetector::Point_2 sc1A(m_InnerRadius, 0), sc2A(0, m_InnerRadius), sc3A(-m_InnerRadius, 0);
  Circle_2 inner_circleA(sc1A, sc2A, sc3A);
  std::vector<CGAL::Object> resA;
  CGAL::intersection(inner_circleA, perpA, std::back_inserter(resA));
  std::vector<CGAL::Object>::const_iterator iterA;
  double pxmax = 0.;
  for (iterA = resA.begin(); iterA != resA.end(); ++iterA)
  {
    CGAL::Object obj = *iterA;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned>>(&obj))
    {
      if (CGAL::to_double(point->first.x()) > pxmax)
      {
        pxmax = CGAL::to_double(point->first.x());
        PHG4OuterHcalDetector::Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
        upperleft = pntmp;
      }
    }
    else
    {
      std::cout << "CGAL::Object type not pair..." << std::endl;
    }
  }
  PHG4OuterHcalDetector::Point_2 so1A(m_OuterRadius, 0), so2A(0, m_OuterRadius), so3A(-m_OuterRadius, 0);
  Circle_2 outer_circleA(so1A, so2A, so3A);
  resA.clear();  // just clear the content from the last intersection search
  CGAL::intersection(outer_circleA, perpA, std::back_inserter(resA));
  for (iterA = resA.begin(); iterA != resA.end(); ++iterA)
  {
    CGAL::Object obj = *iterA;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned>>(&obj))
    {
      if (CGAL::to_double(point->first.x()) > CGAL::to_double(p_loweredge.x()))
      {
        PHG4OuterHcalDetector::Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
        upperright = pntmp;
      }
    }
    else
    {
      std::cout << "CGAL::Object type not pair..." << std::endl;
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
  G4VSolid *steel_plate_uncut = new G4ExtrudedSolid("SteelPlateUnCut",
                                                    vertexes,
                                                    m_SizeZ / 2.0,
                                                    zero, 1.0,
                                                    zero, 1.0);

  m_VolumeSteel = steel_plate_uncut->GetCubicVolume() * m_NumScintiPlates;
  // now cut out space for magnet at the ends
  if (!m_SteelCutoutForMagnetG4Solid)
  {
    return steel_plate_uncut;
  }
  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateX(-90 * deg);
  G4VSolid *steel_firstcut_solid = new G4SubtractionSolid("SteelPlateFirstCut", steel_plate_uncut, m_SteelCutoutForMagnetG4Solid, rotm, G4ThreeVector(0, 0, 0));
  delete rotm;
  //   DisplayVolume(steel_plate_uncut, hcalenvelope);
  //    DisplayVolume(m_SteelCutoutForMagnetG4Solid, hcalenvelope);
  //    DisplayVolume(m_SteelCutoutForMagnetG4Solid, hcalenvelope,rotm);
  //    DisplayVolume(steel_firstcut_solid, hcalenvelope);
  rotm = new G4RotationMatrix();
  rotm->rotateX(90 * deg);
  G4VSolid *steel_cut_solid = new G4SubtractionSolid("SteelPlateCut", steel_firstcut_solid, m_SteelCutoutForMagnetG4Solid, rotm, G4ThreeVector(0, 0, 0));
  //           DisplayVolume(steel_cut_solid, hcalenvelope);
  delete rotm;
  return steel_cut_solid;
}

void PHG4OuterHcalDetector::ShiftSecantToTangent(PHG4OuterHcalDetector::Point_2 &lowleft, PHG4OuterHcalDetector::Point_2 &upleft, PHG4OuterHcalDetector::Point_2 &upright, PHG4OuterHcalDetector::Point_2 &lowright)
{
  Line_2 secant(lowleft, upleft);
  Segment_2 upedge(upleft, upright);
  Segment_2 lowedge(lowleft, lowright);
  double xmid = (CGAL::to_double(lowleft.x()) + CGAL::to_double(upleft.x())) / 2.;
  double ymid = (CGAL::to_double(lowleft.y()) + CGAL::to_double(upleft.y())) / 2.;
  PHG4OuterHcalDetector::Point_2 midpoint(xmid, ymid);
  Line_2 sekperp = secant.perpendicular(midpoint);
  PHG4OuterHcalDetector::Point_2 sc1(m_InnerRadius, 0), sc2(0, m_InnerRadius), sc3(-m_InnerRadius, 0);
  Circle_2 inner_circle(sc1, sc2, sc3);
  std::vector<CGAL::Object> res;
  CGAL::intersection(inner_circle, sekperp, std::back_inserter(res));
  std::vector<CGAL::Object>::const_iterator iter;
  double pxmax = 0.;
  PHG4OuterHcalDetector::Point_2 tangtouch;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned>>(&obj))
    {
      if (CGAL::to_double(point->first.x()) > pxmax)
      {
        pxmax = CGAL::to_double(point->first.x());
        PHG4OuterHcalDetector::Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
        tangtouch = pntmp;
      }
    }
    else
    {
      std::cout << "CGAL::Object type not pair..." << std::endl;
    }
  }
  Line_2 leftside = sekperp.perpendicular(tangtouch);
  CGAL::Object result = CGAL::intersection(upedge, leftside);
  if (const PHG4OuterHcalDetector::Point_2 *ipoint = CGAL::object_cast<PHG4OuterHcalDetector::Point_2>(&result))
  {
    upleft = *ipoint;
  }
  result = CGAL::intersection(lowedge, leftside);
  if (const PHG4OuterHcalDetector::Point_2 *ipoint = CGAL::object_cast<PHG4OuterHcalDetector::Point_2>(&result))
  {
    lowleft = *ipoint;
  }
  return;
}

void PHG4OuterHcalDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
#ifdef SCINTITEST
  ConstructOuterHcal(logicWorld);
  return;
#endif
  recoConsts *rc = recoConsts::instance();
  G4Material *Air = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  G4VSolid *hcal_envelope_cylinder = new G4Tubs("OuterHcal_envelope_solid", m_EnvelopeInnerRadius, m_EnvelopeOuterRadius, m_EnvelopeZ / 2., 0, 2 * M_PI);
  m_VolumeEnvelope = hcal_envelope_cylinder->GetCubicVolume();
  G4LogicalVolume *hcal_envelope_log = new G4LogicalVolume(hcal_envelope_cylinder, Air, G4String("OuterHcal_envelope"), nullptr, nullptr, nullptr);
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(m_Params->get_double_param("rot_x") * deg);
  hcal_rotm.rotateY(m_Params->get_double_param("rot_y") * deg);
  hcal_rotm.rotateZ(m_Params->get_double_param("rot_z") * deg);
  G4VPhysicalVolume *mothervol = new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm)), hcal_envelope_log, "OuterHcal", logicWorld, false, false, OverlapCheck());
  m_DisplayAction->SetMyTopVolume(mothervol);
  ConstructOuterHcal(hcal_envelope_log);
  std::vector<G4VPhysicalVolume *>::iterator it = m_ScintiMotherAssembly->GetVolumesIterator();
  for (unsigned int i = 0; i < m_ScintiMotherAssembly->TotalImprintedVolumes(); i++)
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
    // use boost tokenizer to separate the _, then take value
    // after "impr" for mother volume and after "pv" for scintillator slat
    // use boost lexical cast for string -> int conversion
    // the CopyNo is the mother volume + scinti id
    // so we can use the CopyNo rather than decoding the string further
    // looking for "pv"
    boost::char_separator<char> sep("_");
    boost::tokenizer<boost::char_separator<char>> tok((*it)->GetName(), sep);
    boost::tokenizer<boost::char_separator<char>>::const_iterator tokeniter;
    for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
    {
      if (*tokeniter == "impr")
      {
        ++tokeniter;
        if (tokeniter != tok.end())
        {
          int layer_id = boost::lexical_cast<int>(*tokeniter);
          // check detector description, for assemblyvolumes it is not possible
          // to give the first volume id=0, so they go from id=1 to id=n.
          // I am not going to start with fortran again - our indices start
          // at zero, id=0 to id=n-1. So subtract one here
          int tower_id = (*it)->GetCopyNo() - layer_id;
          layer_id--;
          std::pair<int, int> layer_twr = std::make_pair(layer_id, tower_id);
          m_ScintiTilePhysVolMap.insert(std::pair<G4VPhysicalVolume *, std::pair<int, int>>(*it, layer_twr));
          if (layer_id < 0 || layer_id >= m_NumScintiPlates)
          {
            std::cout << "invalid scintillator row " << layer_id
                      << ", valid range 0 < row < " << m_NumScintiPlates << std::endl;
            gSystem->Exit(1);
          }
        }
        else
        {
          std::cout << PHWHERE << " Error parsing " << (*it)->GetName()
                    << " for mother volume number " << std::endl;
          gSystem->Exit(1);
        }
        break;
      }
    }
    ++it;
  }
  return;
}

int PHG4OuterHcalDetector::ConstructOuterHcal(G4LogicalVolume *hcalenvelope)
{
  ConsistencyCheck();
  SetTiltViaNcross();  // if number of crossings is set, use it to determine tilt
  CheckTiltAngle();    // die if the tilt angle is out of range
  // the needed steel cutout volume for the magnet is constructed with
  // the scintillators since we have the theta angle at that point

  // call field setup here where we have the calculated tilt angle if number
  // of crossings is given
  m_FieldSetup = new PHG4OuterHcalFieldSetup(
      m_NumScintiPlates, /*G4int steelPlates*/
      m_ScintiGap,       /*G4double scintiGap*/
      m_TiltAngle);      /*G4double tiltAngle*/

  m_ScintiMotherAssembly = ConstructHcalScintillatorAssembly(hcalenvelope);
#ifdef SCINTITEST
  return 0;
#endif
  G4VSolid *steel_plate = ConstructSteelPlate(hcalenvelope);
  //   DisplayVolume(steel_plate_4 ,hcalenvelope);
  G4LogicalVolume *steel_logical = new G4LogicalVolume(steel_plate, GetDetectorMaterial(m_Params->get_string_param("material")), "HcalOuterSteelPlate", nullptr, nullptr, nullptr);
  m_DisplayAction->AddSteelVolume(steel_logical);
  double phi = 0;
  double deltaphi = 2 * M_PI / m_NumScintiPlates;
  std::ostringstream name;
  double middlerad = m_OuterRadius - (m_OuterRadius - m_InnerRadius) / 2.;
  // okay this is crude. Since the inner and outer radius of the scintillator is different from the inner/outer
  // radius of the steel so the scintillator needs some shifting to get the gaps right
  // basically the first shift (sumshift) puts it at the inner radius of the steel
  // then it needs to be shifted up by twice the difference between the inner steel and inner scintillator radius
  // since sumshift is div by 2 for G4 but we need to shift the full delta(inner radius)
  double scinti_tile_orig_length = m_ScintiTileXUpper + m_ScintiTileXLower - subtract_from_scinti_x;
  double shiftup = (m_ScintiInnerRadius - m_InnerRadius) / cos(m_TiltAngle / rad);
  double sumshift = scinti_tile_orig_length - m_ScintiTileX;
  sumshift = sumshift - 2 * shiftup;
  double shiftslat = fabs(m_ScintiTileXLower - m_ScintiTileXUpper) / 2. + sumshift / 2.;
  // calculate phi offset (copied from code inside following loop):
  // first get the center point (phi=0) so it's middlerad/0
  // then shift the scintillator center as documented in loop
  // then
  // for positive tilt angles we need the lower left corner of the scintillator
  // for negative tilt angles we nee the upper right corner of the scintillator
  // as it turns out the code uses the middle of the face of the scintillator
  // as reference, if this is a problem the code needs to be modified to
  // actually calculate the corner (but the math of the construction is that
  // the middle of the scintillator sits at zero)
  double xp = cos(phi) * middlerad;
  double yp = sin(phi) * middlerad;
  xp -= cos((-m_TiltAngle) / rad - phi) * shiftslat;
  yp += sin((-m_TiltAngle) / rad - phi) * shiftslat;
  if (m_TiltAngle > 0)
  {
    double xo = xp - (m_ScintiTileX / 2.) * cos(m_TiltAngle / rad);
    double yo = yp - (m_ScintiTileX / 2.) * sin(m_TiltAngle / rad);
    phi = -atan(yo / xo);
  }
  else if (m_TiltAngle < 0)
  {
    double xo = xp + (m_ScintiTileX / 2.) * cos(m_TiltAngle / rad);
    double yo = yp + (m_ScintiTileX / 2.) * sin(m_TiltAngle / rad);
    phi = -atan(yo / xo);
  }
  // else (for m_TiltAngle = 0) phi stays zero
  for (int i = 0; i < m_NumScintiPlates; i++)
  {
    G4RotationMatrix *Rot = new G4RotationMatrix();
    double ypos = sin(phi) * middlerad;
    double xpos = cos(phi) * middlerad;
    // the center of the scintillator is not the center of the inner hcal
    // but depends on the tilt angle. Therefore we need to shift
    // the center from the mid point
    ypos += sin((-m_TiltAngle) / rad - phi) * shiftslat;
    xpos -= cos((-m_TiltAngle) / rad - phi) * shiftslat;
    Rot->rotateZ(phi * rad + m_TiltAngle);
    G4ThreeVector g4vec(xpos, ypos, 0);
    // great the MakeImprint always adds 1 to the copy number and 0 has a
    // special meaning (which then also adds 1). Basically our volume names
    // will start at 1 instead of 0 and there is nothing short of patching this
    // method. I'll take care of this in the decoding of the volume name
    // AAAAAAARHGS
    m_ScintiMotherAssembly->MakeImprint(hcalenvelope, g4vec, Rot, i, OverlapCheck());
    delete Rot;
    Rot = new G4RotationMatrix();
    Rot->rotateZ(-phi * rad);
    name.str("");
    name << "OuterHcalSteel_" << i;
    m_SteelAbsorberVec.insert(new G4PVPlacement(Rot, G4ThreeVector(0, 0, 0), steel_logical, name.str(), hcalenvelope, false, i, OverlapCheck()));
    phi += deltaphi;
  }
  hcalenvelope->SetFieldManager(m_FieldSetup->get_Field_Manager_Gap(), false);

  steel_logical->SetFieldManager(m_FieldSetup->get_Field_Manager_Iron(), true);
  return 0;
}

void PHG4OuterHcalDetector::ConstructHcalSingleScintillators(G4LogicalVolume *hcalenvelope)
{
  G4VSolid *bigtile = ConstructScintillatorBox(hcalenvelope);
  // eta->theta
  G4double delta_eta = m_Params->get_double_param("scinti_eta_coverage") / m_NumScintiTiles;
  G4double eta = 0;
  G4double theta;
  G4double x[4];
  G4double z[4];
  std::ostringstream name;
  double overhang = (m_ScintiTileX - (m_ScintiOuterRadius - m_ScintiInnerRadius)) / 2.;
  double offset = 1 * cm + overhang;  // add 1cm to make sure the G4ExtrudedSolid
  // is larger than the tile so we do not have
  // funny edge effects when overlapping vols
  double magnet_cutout_x = (m_Params->get_double_param("magnet_cutout_scinti_radius") * cm - m_ScintiInnerRadius) / cos(m_TiltAngle / rad);
  double x_inner = m_ScintiInnerRadius - overhang;
  double inner_offset = offset;
  // coordinates like the steel plates:
  // 0/0 upper left
  // 1/1 upper right
  // 2/2 lower right
  // 3/3 lower left
  // sorry they are different than the coordinates used for the scintilators
  // here, this is why the indices are seemingly mixed up
  double xsteelcut[4];
  double zsteelcut[4];
  std::fill_n(zsteelcut, 4, NAN);
  double steel_overhang = (m_ScintiTileXUpper + m_ScintiTileXLower - subtract_from_scinti_x - (m_OuterRadius - m_InnerRadius)) / 2.;
  double steel_offset = 1 * cm + steel_overhang;  // add 1cm to make sure the G4ExtrudedSolid
  double steel_x_inner = m_InnerRadius - steel_overhang;
  double steel_magnet_cutout_x = (m_Params->get_double_param("magnet_cutout_radius") * cm - m_InnerRadius) / cos(m_TiltAngle / rad);
  double steel_inner_offset = steel_offset;
  xsteelcut[0] = steel_x_inner + steel_magnet_cutout_x;
  xsteelcut[1] = xsteelcut[0];
  xsteelcut[2] = m_InnerRadius - steel_offset;
  xsteelcut[3] = xsteelcut[2];
  double scinti_gap_neighbor = m_Params->get_double_param("scinti_gap_neighbor") * cm;
  for (int i = 0; i < m_NumScintiTiles; i++)
  {
    if (i >= m_Params->get_int_param("magnet_cutout_first_scinti"))
    {
      x_inner = m_ScintiInnerRadius - overhang + magnet_cutout_x;
      inner_offset = offset - magnet_cutout_x;
    }
    theta = M_PI / 2 - PHG4Utils::get_theta(eta);  // theta = 90 for eta=0
    x[0] = x_inner;
    z[0] = tan(theta) * m_ScintiInnerRadius;
    x[1] = m_ScintiOuterRadius + overhang;  // since the tile is tilted, x is not at the outer radius but beyond
    z[1] = tan(theta) * m_ScintiOuterRadius;
    if (i >= m_Params->get_int_param("magnet_cutout_first_scinti"))
    {
      z[0] = tan(theta) * (m_ScintiInnerRadius + (m_Params->get_double_param("magnet_cutout_scinti_radius") * cm - m_ScintiInnerRadius));
    }
    eta += delta_eta;
    theta = M_PI / 2 - PHG4Utils::get_theta(eta);  // theta = 90 for eta=0
    x[2] = x_inner;
    z[2] = tan(theta) * m_ScintiInnerRadius;
    if (i >= m_Params->get_int_param("magnet_cutout_first_scinti"))
    {
      z[2] = tan(theta) * (m_ScintiInnerRadius + (m_Params->get_double_param("magnet_cutout_scinti_radius") * cm - m_ScintiInnerRadius));
    }
    x[3] = m_ScintiOuterRadius + overhang;  // since the tile is tilted, x is not at the outer radius but beyond
    z[3] = tan(theta) * m_ScintiOuterRadius;
    // apply gap between scintillators
    z[0] += scinti_gap_neighbor / 2.;
    z[1] += scinti_gap_neighbor / 2.;
    z[2] -= scinti_gap_neighbor / 2.;
    z[3] -= scinti_gap_neighbor / 2.;
    PHG4OuterHcalDetector::Point_2 leftsidelow(z[0], x[0]);
    PHG4OuterHcalDetector::Point_2 leftsidehigh(z[1], x[1]);
    x[0] = m_ScintiInnerRadius - inner_offset;
    z[0] = x_at_y(leftsidelow, leftsidehigh, x[0]);
    x[1] = m_ScintiOuterRadius + offset;
    z[1] = x_at_y(leftsidelow, leftsidehigh, x[1]);
    PHG4OuterHcalDetector::Point_2 rightsidelow(z[2], x[2]);
    PHG4OuterHcalDetector::Point_2 rightsidehigh(z[3], x[3]);
    x[2] = m_ScintiOuterRadius + offset;
    z[2] = x_at_y(rightsidelow, rightsidehigh, x[2]);
    x[3] = m_ScintiInnerRadius - inner_offset;
    z[3] = x_at_y(rightsidelow, rightsidehigh, x[3]);
    // store corner points of extruded solid we need to subtract from steel
    if (i == m_Params->get_int_param("magnet_cutout_first_scinti"))
    {
      double x0 = m_InnerRadius - (steel_inner_offset - steel_magnet_cutout_x);
      double z0 = x_at_y(leftsidelow, leftsidehigh, x0);
      double xpos = m_InnerRadius - steel_offset;
      zsteelcut[0] = z0;
      zsteelcut[3] = x_at_y(leftsidelow, leftsidehigh, xpos);
    }
    double x2 = m_OuterRadius + steel_offset;
    double z2 = x_at_y(rightsidelow, rightsidehigh, x2);
    zsteelcut[1] = z2 + 1 * cm;
    zsteelcut[2] = z2 + 1 * cm;
    std::vector<G4TwoVector> vertexes;
    for (int j = 0; j < 4; j++)
    {
      G4TwoVector v(x[j], z[j]);
      vertexes.push_back(v);
    }
    G4TwoVector zero(0, 0);

    G4VSolid *scinti = new G4ExtrudedSolid("ScintillatorTile",
                                           vertexes,
                                           m_ScintiTileThickness + 0.2 * mm,
                                           zero, 1.0,
                                           zero, 1.0);
    G4RotationMatrix *rotm = new G4RotationMatrix();
    rotm->rotateX(-90 * deg);
    name.str("");
    name << "scintillator_" << i << "_left";
    G4VSolid *scinti_tile = new G4IntersectionSolid(name.str(), bigtile, scinti, rotm, G4ThreeVector(-(m_ScintiInnerRadius + m_ScintiOuterRadius) / 2., 0, 0));
    delete rotm;
    m_ScintiTilesVec[i + m_NumScintiTiles] = scinti_tile;
    rotm = new G4RotationMatrix();
    rotm->rotateX(90 * deg);
    name.str("");
    name << "scintillator_" << i << "_right";
    scinti_tile = new G4IntersectionSolid(name.str(), bigtile, scinti, rotm, G4ThreeVector(-(m_ScintiInnerRadius + m_ScintiOuterRadius) / 2., 0, 0));
    delete rotm;
    m_ScintiTilesVec[m_NumScintiTiles - i - 1] = scinti_tile;
  }
#ifdef SCINTITEST
  for (unsigned int i = 0; i < m_ScintiTilesVec.size(); i++)
  {
    if (m_ScintiTilesVec[i])
    {
      DisplayVolume(m_ScintiTilesVec[i], hcalenvelope);
    }
  }
#endif

  std::vector<G4TwoVector> vertexes;
  for (int j = 0; j < 4; j++)
  {
    if (!isfinite(zsteelcut[j]))
    {
      return;
    }
    G4TwoVector v(xsteelcut[j], zsteelcut[j]);
    vertexes.push_back(v);
  }
  G4TwoVector zero(0, 0);
  m_SteelCutoutForMagnetG4Solid = new G4ExtrudedSolid("ScintillatorTile",
                                                      vertexes,
                                                      m_ScintiTileThickness + 20 * cm,
                                                      zero, 1.0,
                                                      zero, 1.0);
  return;
}

G4double
PHG4OuterHcalDetector::x_at_y(PHG4OuterHcalDetector::Point_2 &p0, PHG4OuterHcalDetector::Point_2 &p1, G4double yin)
{
  double xret = NAN;
  double x[2];
  x[0] = CGAL::to_double(p0.x());
  x[1] = CGAL::to_double(p1.x());
  Line_2 l(p0, p1);
  double newx = fabs(x[0]) + fabs(x[1]);
  PHG4OuterHcalDetector::Point_2 p0new(-newx, yin);
  PHG4OuterHcalDetector::Point_2 p1new(newx, yin);
  Segment_2 seg(p0new, p1new);
  CGAL::Object result = CGAL::intersection(l, seg);
  if (const PHG4OuterHcalDetector::Point_2 *ipoint = CGAL::object_cast<PHG4OuterHcalDetector::Point_2>(&result))
  {
    xret = CGAL::to_double(ipoint->x());
  }
  else
  {
    std::cout << PHWHERE << " failed for y = " << yin << std::endl;
    std::cout << "p0(x): " << CGAL::to_double(p0.x()) << ", p0(y): " << CGAL::to_double(p0.y()) << std::endl;
    std::cout << "p1(x): " << CGAL::to_double(p1.x()) << ", p1(y): " << CGAL::to_double(p1.y()) << std::endl;
    exit(1);
  }
  return xret;
}

G4AssemblyVolume *
PHG4OuterHcalDetector::ConstructHcalScintillatorAssembly(G4LogicalVolume *hcalenvelope)
{
#ifdef SCINTITEST
  ConstructHcalSingleScintillators(hcalenvelope);
  return nullptr;
#endif
  ConstructHcalSingleScintillators(hcalenvelope);
  G4AssemblyVolume *assmeblyvol = new G4AssemblyVolume();
  std::ostringstream name;
  G4ThreeVector g4vec;
  double steplimits = m_Params->get_double_param("steplimits") * cm;
  for (unsigned int i = 0; i < m_ScintiTilesVec.size(); i++)
  {
    name.str("");
    name << m_ScintiLogicNamePrefix << i;
    G4UserLimits *g4userlimits = nullptr;
    if (isfinite(steplimits))
    {
      g4userlimits = new G4UserLimits(steplimits);
    }
    G4LogicalVolume *scinti_tile_logic = new G4LogicalVolume(m_ScintiTilesVec[i], GetDetectorMaterial("G4_POLYSTYRENE"), name.str(), nullptr, nullptr, g4userlimits);
    m_DisplayAction->AddScintiVolume(scinti_tile_logic);
    assmeblyvol->AddPlacedVolume(scinti_tile_logic, g4vec, nullptr);

    //field after burner
    scinti_tile_logic->SetFieldManager(m_FieldSetup->get_Field_Manager_Gap(), true);
  }
  return assmeblyvol;
}

int PHG4OuterHcalDetector::ConsistencyCheck() const
{
  // just make sure the parameters make a bit of sense
  if (m_InnerRadius >= m_OuterRadius)
  {
    std::cout << PHWHERE << ": Inner Radius " << m_InnerRadius / cm
              << " cm larger than Outer Radius " << m_OuterRadius / cm
              << " cm" << std::endl;
    gSystem->Exit(1);
  }
  if (m_ScintiTileThickness > m_ScintiGap)
  {
    std::cout << PHWHERE << "Scintillator thickness " << m_ScintiTileThickness / cm
              << " cm larger than scintillator gap " << m_ScintiGap / cm
              << " cm" << std::endl;
    gSystem->Exit(1);
  }
  if (m_ScintiOuterRadius <= m_ScintiInnerRadius)
  {
    std::cout << PHWHERE << "Scintillator outer radius " << m_ScintiOuterRadius / cm
              << " cm smaller than scintillator inner radius " << m_ScintiInnerRadius / cm
              << " cm" << std::endl;
    gSystem->Exit(1);
  }
  if (m_ScintiOuterRadius <= m_InnerRadius)
  {
    std::cout << PHWHERE << "Scintillator outer radius " << m_ScintiOuterRadius / cm
              << " cm smaller than inner radius " << m_InnerRadius / cm
              << " cm" << std::endl;
    gSystem->Exit(1);
  }
  if (m_ScintiInnerRadius >= m_OuterRadius)
  {
    std::cout << PHWHERE << "Scintillator inner radius " << m_ScintiInnerRadius / cm
              << " cm larger than inner radius " << m_InnerRadius / cm
              << " cm" << std::endl;
    gSystem->Exit(1);
  }
  if (m_Params->get_double_param("magnet_cutout_scinti_radius") * cm < m_ScintiInnerRadius)
  {
    std::cout << PHWHERE << "Magnet scintillator cutout radius " << m_Params->get_double_param("magnet_cutout_scinti_radius")
              << " cm smaller than inner scintillator radius " << m_ScintiInnerRadius / cm
              << " cm" << std::endl;
    gSystem->Exit(1);
  }
  if (m_Params->get_double_param("magnet_cutout_radius") * cm < m_InnerRadius)
  {
    std::cout << PHWHERE << "Magnet steel cutout radius " << m_Params->get_double_param("magnet_cutout_radius")
              << " cm smaller than inner radius " << m_InnerRadius / cm
              << " cm" << std::endl;
    gSystem->Exit(1);
  }

  return 0;
}

void PHG4OuterHcalDetector::SetTiltViaNcross()
{
  int ncross = m_Params->get_int_param("ncross");
  if (!ncross || isfinite(m_TiltAngle))
  {
    // mark ncross parameter as not used
    m_Params->set_int_param("ncross", 0);
    return;
  }
  if ((isfinite(m_TiltAngle)) && (Verbosity() > 0))
  {
    std::cout << "both number of crossings and tilt angle are set" << std::endl;
    std::cout << "using number of crossings to determine tilt angle" << std::endl;
  }
  double mid_radius = m_InnerRadius + (m_OuterRadius - m_InnerRadius) / 2.;
  double deltaphi = (2 * M_PI / m_NumScintiPlates) * ncross;
  PHG4OuterHcalDetector::Point_2 pnull(0, 0);
  PHG4OuterHcalDetector::Point_2 plow(m_InnerRadius, 0);
  PHG4OuterHcalDetector::Point_2 phightmp(1, tan(deltaphi));
  PHG4OuterHcalDetector::Point_2 pin1(m_InnerRadius, 0), pin2(0, m_InnerRadius), pin3(-m_InnerRadius, 0);
  Circle_2 inner_circle(pin1, pin2, pin3);
  PHG4OuterHcalDetector::Point_2 pmid1(mid_radius, 0), pmid2(0, mid_radius), pmid3(-mid_radius, 0);
  Circle_2 mid_circle(pmid1, pmid2, pmid3);
  PHG4OuterHcalDetector::Point_2 pout1(m_OuterRadius, 0), pout2(0, m_OuterRadius), pout3(-m_OuterRadius, 0);
  Circle_2 outer_circle(pout1, pout2, pout3);
  Line_2 l_up(pnull, phightmp);
  std::vector<CGAL::Object> res;
  CGAL::intersection(outer_circle, l_up, std::back_inserter(res));
  PHG4OuterHcalDetector::Point_2 upperright;
  std::vector<CGAL::Object>::const_iterator iter;
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned>>(&obj))
    {
      if (CGAL::to_double(point->first.x()) > 0)
      {
        PHG4OuterHcalDetector::Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
        upperright = pntmp;
      }
    }
    else
    {
      std::cout << "CGAL::Object type not pair..." << std::endl;
      exit(1);
    }
  }
  Line_2 l_right(plow, upperright);
  res.clear();
  PHG4OuterHcalDetector::Point_2 midpoint;
  CGAL::intersection(mid_circle, l_right, std::back_inserter(res));
  for (iter = res.begin(); iter != res.end(); ++iter)
  {
    CGAL::Object obj = *iter;
    if (const std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned> *point = CGAL::object_cast<std::pair<CGAL::Circular_arc_point_2<PHG4OuterHcalDetector::Circular_k>, unsigned>>(&obj))
    {
      if (CGAL::to_double(point->first.x()) > 0)
      {
        PHG4OuterHcalDetector::Point_2 pntmp(CGAL::to_double(point->first.x()), CGAL::to_double(point->first.y()));
        midpoint = pntmp;
      }
    }
    else
    {
      std::cout << "CGAL::Object type not pair..." << std::endl;
      exit(1);
    }
  }
  // length left side
  double ll = sqrt((CGAL::to_double(midpoint.x()) - m_InnerRadius) * (CGAL::to_double(midpoint.x()) - m_InnerRadius) + CGAL::to_double(midpoint.y()) * CGAL::to_double(midpoint.y()));
  double upside = sqrt(CGAL::to_double(midpoint.x()) * CGAL::to_double(midpoint.x()) + CGAL::to_double(midpoint.y()) * CGAL::to_double(midpoint.y()));
  //  c^2 = a^2+b^2 - 2ab*cos(gamma)
  // gamma = acos((a^2+b^2=c^2)/2ab
  double tiltangle = acos((ll * ll + upside * upside - m_InnerRadius * m_InnerRadius) / (2 * ll * upside));
  tiltangle = tiltangle * rad;
  m_TiltAngle = copysign(tiltangle, ncross);
  m_Params->set_double_param("tilt_angle", m_TiltAngle / deg);
  return;
}

// check if tilt angle is reasonable - too large, no intersections with inner radius
int PHG4OuterHcalDetector::CheckTiltAngle() const
{
  if (fabs(m_TiltAngle) >= M_PI)
  {
    std::cout << PHWHERE << "invalid tilt angle, abs(tilt) >= 90 deg: " << (m_TiltAngle / deg)
              << std::endl;
    exit(1);
  }

  double mid_radius = m_InnerRadius + (m_OuterRadius - m_InnerRadius) / 2.;
  PHG4OuterHcalDetector::Point_2 pmid(mid_radius, 0);  // center of scintillator
  double xcoord = 0;
  double ycoord = mid_radius * tan(m_TiltAngle / rad);
  PHG4OuterHcalDetector::Point_2 pxnull(xcoord, ycoord);
  Line_2 s2(pmid, pxnull);
  PHG4OuterHcalDetector::Point_2 sc1(m_InnerRadius, 0), sc2(0, m_InnerRadius), sc3(-m_InnerRadius, 0);
  Circle_2 inner_circle(sc1, sc2, sc3);
  std::vector<CGAL::Object> res;
  CGAL::intersection(inner_circle, s2, std::back_inserter(res));
  if (res.size() == 0)
  {
    std::cout << PHWHERE << " Tilt angle " << (m_TiltAngle / deg)
              << " too large, no intersection with inner radius" << std::endl;
    exit(1);
  }
  return 0;
}

void PHG4OuterHcalDetector::Print(const std::string &what) const
{
  std::cout << "Outer Hcal Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Volume Envelope: " << m_VolumeEnvelope / cm / cm / cm << " cm^3" << std::endl;
    std::cout << "Volume Steel: " << m_VolumeSteel / cm / cm / cm << " cm^3" << std::endl;
    std::cout << "Volume Scintillator: " << m_VolumeScintillator / cm / cm / cm << " cm^3" << std::endl;
    std::cout << "Volume Air: " << (m_VolumeEnvelope - m_VolumeSteel - m_VolumeScintillator) / cm / cm / cm << " cm^3" << std::endl;
  }
  return;
}

std::pair<int, int> PHG4OuterHcalDetector::GetLayerTowerId(G4VPhysicalVolume *volume) const
{
  auto it = m_ScintiTilePhysVolMap.find(volume);
  if (it != m_ScintiTilePhysVolMap.end())
  {
    return it->second;
  }
  std::cout << "could not locate volume " << volume->GetName()
            << " in Inner Hcal scintillator map" << std::endl;
  gSystem->Exit(1);
  // that's dumb but code checkers do not know that gSystem->Exit()
  // terminates, so using the standard exit() makes them happy
  exit(1);
}
