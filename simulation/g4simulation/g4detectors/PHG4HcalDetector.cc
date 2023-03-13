#include "PHG4HcalDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeomv3.h"

#include <g4main/PHG4Detector.h>  // for PHG4Detector
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Geant4/G4Colour.hh>  // for G4Colour
#include <Geant4/G4Cons.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>  // for cm, deg
#include <Geant4/G4ThreeVector.hh>    // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>    // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4VisAttributes.hh>

#include <algorithm>  // for max
#include <cmath>      // for sin, cos, sqrt, M_PI, asin
#include <cstdlib>    // for exit
#include <iostream>   // for operator<<, basic_ostream
#include <sstream>
#include <utility>  // for pair
#include <vector>   // for vector

class PHG4CylinderGeom;

// uncomment if you want to make a graphics display where the slats are visible
// it makes them stick out of the hcal for visibility
// NEVER EVER RUN REAL SIMS WITH THIS
//#define DISPLAY

int PHG4HcalDetector::INACTIVE = -100;
//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4HcalDetector::PHG4HcalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , TrackerMaterial(nullptr)
  , TrackerThickness(100 * cm)
  , cylinder_logic(nullptr)
  , cylinder_physi(nullptr)
  , radius(100 * cm)
  , length(100 * cm)
  , xpos(0 * cm)
  , ypos(0 * cm)
  , zpos(0 * cm)
  , _sciTilt(0)
  , _sciWidth(0.6 * cm)
  , _sciNum(100)
  , _sciPhi0(0)
  , _region(nullptr)
  , active(0)
  , absorberactive(0)
  , layer(lyr)
{
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4HcalDetector::IsInCylinderActive(const G4VPhysicalVolume* volume)
{
  //  std::cout << "checking detector" << std::endl;
  if (active && box_vol.find(volume) != box_vol.end())
  {
    return box_vol.find(volume)->second;
  }
  if (absorberactive && volume == cylinder_physi)
  {
    return -1;
  }
  return INACTIVE;
}

//_______________________________________________________________
void PHG4HcalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  TrackerMaterial = GetDetectorMaterial(material);

  G4Tubs* _cylinder_solid = new G4Tubs(G4String(GetName()),
                                       radius,
                                       radius + TrackerThickness,
                                       length / 2.0, 0, twopi);
  double innerlength = PHG4Utils::GetLengthForRapidityCoverage(radius) * 2;
  double deltalen = (length - innerlength) / 2.;  // length difference on one side
  double cone_size_multiplier = 1.01;             // 1 % larger
  double cone_thickness = TrackerThickness * cone_size_multiplier;
  double inner_cone_radius = radius - ((cone_thickness - TrackerThickness) / 2.);
  double cone_length = deltalen * cone_size_multiplier;
  G4Cons* cone2 = new G4Cons("conehead2",
                             inner_cone_radius, inner_cone_radius,
                             inner_cone_radius, inner_cone_radius + cone_thickness,
                             cone_length / 2.0, 0, twopi);
  G4Cons* cone1 = new G4Cons("conehead",
                             inner_cone_radius, inner_cone_radius + cone_thickness,
                             inner_cone_radius, inner_cone_radius,
                             cone_length / 2.0, 0, twopi);
  double delta_len = cone_length - deltalen;
  G4ThreeVector zTransneg(0, 0, -(length - cone_length + delta_len) / 2.0);
  G4ThreeVector zTranspos(0, 0, (length - cone_length + delta_len) / 2.0);
  G4SubtractionSolid* subtraction_tmp =
      new G4SubtractionSolid("Cylinder-Cone", _cylinder_solid, cone1, nullptr, zTransneg);
  G4SubtractionSolid* subtraction =
      new G4SubtractionSolid("Cylinder-Cone-Cone", subtraction_tmp, cone2, nullptr, zTranspos);
  cylinder_logic = new G4LogicalVolume(subtraction,
                                       TrackerMaterial,
                                       G4String(GetName()),
                                       nullptr, nullptr, nullptr);
  G4VisAttributes* VisAtt = new G4VisAttributes();
  VisAtt->SetColour(G4Colour::Grey());
  VisAtt->SetVisibility(true);
  VisAtt->SetForceSolid(true);
  cylinder_logic->SetVisAttributes(VisAtt);

  cylinder_physi = new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, zpos),
                                     cylinder_logic,
                                     G4String(GetName()),
                                     logicWorld, false, false, OverlapCheck());
  // Figure out corners of scintillator inside the containing G4Tubs.
  // Work our way around the scintillator cross section in a counter
  // clockwise fashion: ABCD
  double r1 = radius;
  double r2 = radius + TrackerThickness;

  // The coordinates of the inner corners of the scintillator
  double x4 = r1;
  double y4 = _sciWidth / 2.0;

  double x1 = r1;
  double y1 = -y4;

  double a = _sciTilt * M_PI / 180.0;

  // The parametric equation for the line from A along the side of the
  // scintillator is (x,y) = (x1,y1) + u * (cos(a), sin(a))
  double A = 1.0;
  double B = 2 * (x1 * cos(a) + y1 * sin(a));
  double C = x1 * x1 + y1 * y1 - r2 * r2;
  double D = B * B - 4 * A * C;

  // The only sensible solution, given our definitions, is u > 0.
  double u = (-B + sqrt(D)) / 2 * A;

  // Now we can determine one of the outer corners
  double x2 = x1 + u * cos(a);
  double y2 = y1 + u * sin(a);

  // Similar procedure for (x3,y3) as for (x2,y2)
  A = 1.0;
  B = 2 * (x4 * cos(a) + y4 * sin(a));
  C = x4 * x4 + y4 * y4 - r2 * r2;
  D = B * B - 4 * A * C;
  u = (-B + sqrt(D)) / 2 * A;

  double x3 = x4 + u * cos(a);
  double y3 = y4 + u * sin(a);

  // Now we've got a four-sided "z-section".
  G4TwoVector v1(x1, y1);
  G4TwoVector v2(x2, y2);
  G4TwoVector v3(x3, y3);
  G4TwoVector v4(x4, y4);

  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(v1);
  vertexes.push_back(v2);
  vertexes.push_back(v3);
  vertexes.push_back(v4);

  G4TwoVector zero(0, 0);
  // if you want to make displays where the structure of the hcal is visible
  // add 20 cm to the length of the scintillators
#ifdef DISPLAY
  double blength = length + 20;
#else
  double blength = length;
#endif
  G4ExtrudedSolid* _box_solid = new G4ExtrudedSolid("_BOX",
                                                    vertexes,
                                                    blength / 2.0,
                                                    zero, 1.0,
                                                    zero, 1.0);

  //  double boxlen_half = GetLength(_sciTilt * M_PI / 180.);
  //  G4Box* _box_solid = new G4Box("_BOX", boxlen_half, _sciWidth / 2.0, length / 2.0);

  G4Material* boxmat = GetDetectorMaterial("G4_POLYSTYRENE");
  G4SubtractionSolid* subtractionbox_tmp =
      new G4SubtractionSolid("Box-Cone", _box_solid, cone1, nullptr, zTransneg);
  G4SubtractionSolid* subtractionbox =
      new G4SubtractionSolid("Box-Cone-Cone", subtractionbox_tmp, cone2, nullptr, zTranspos);
  G4LogicalVolume* box_logic = new G4LogicalVolume(subtractionbox,
                                                   boxmat, G4String("BOX"),
                                                   nullptr, nullptr, nullptr);
  VisAtt = new G4VisAttributes();
  PHG4Utils::SetColour(VisAtt, "G4_POLYSTYRENE");
  VisAtt->SetVisibility(true);
  VisAtt->SetForceSolid(true);
  box_logic->SetVisAttributes(VisAtt);

  double phi_increment = 360. / _sciNum;
  std::ostringstream slatname;
  for (int i = 0; i < _sciNum; i++)
  {
    double phi = (i + _sciPhi0) * phi_increment;
    G4ThreeVector myTrans = G4ThreeVector(0, 0, 0);
    G4RotationMatrix Rot(0, 0, 0);
    Rot.rotateZ(phi * deg);
    slatname.str("");
    slatname << "SLAT_" << i;
    G4VPhysicalVolume* box_vol_tmp = new G4PVPlacement(G4Transform3D(Rot, G4ThreeVector(myTrans)),
                                                       box_logic,
                                                       G4String(slatname.str()),
                                                       cylinder_logic, false, false, OverlapCheck());
    box_vol[box_vol_tmp] = i;
  }
  if (active)
  {
    std::ostringstream geonode;
    if (superdetector != "NONE")
    {
      geonode << "CYLINDERGEOM_" << superdetector;
    }
    else
    {
      geonode << "CYLINDERGEOM_" << detector_type << "_" << layer;
    }
    PHG4CylinderGeomContainer* geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode.str());
    if (!geo)
    {
      geo = new PHG4CylinderGeomContainer();
      PHNodeIterator iter(topNode());
      PHCompositeNode* runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
      PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(geo, geonode.str(), "PHObject");
      runNode->addNode(newNode);
    }
    // here in the detector class we have internal units, convert to cm
    // before putting into the geom object
    PHG4CylinderGeom* mygeom = new PHG4CylinderGeomv3(radius / cm, (zpos - length / 2.) / cm, (zpos + length / 2.) / cm, TrackerThickness / cm, _sciNum, _sciTilt * M_PI / 180.0, _sciPhi0 * M_PI / 180.0);
    geo->AddLayerGeom(layer, mygeom);
    //    geo->identify();
  }
}

double
PHG4HcalDetector::GetLength(const double phi) const
{
  double c = radius + TrackerThickness / 2.;
  double b = radius;
  double singamma = sin(phi) * c / b;
  double gamma = M_PI - asin(singamma);
  double alpha = M_PI - gamma - phi;
  double a = c * sin(alpha) / singamma;
  return a;
}

void PHG4HcalDetector::Print(const std::string& /*what*/) const
{
  std::cout << "radius: " << radius << std::endl;
  return;
}
