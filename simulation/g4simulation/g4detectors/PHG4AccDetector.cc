#include "PHG4AccDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeomv2.h"

#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4Utils.h>

#include <fun4all/getClass.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4EllipticalTube.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <sstream>

using namespace std;

int PHG4AccDetector::INACTIVE = -100;
//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4AccDetector::PHG4AccDetector( PHCompositeNode *Node, const std::string &dnam, const int lyr ):
  PHG4Detector(Node, dnam),
  TrackerThickness(100*cm),
  radius(100*cm),
  length(100*cm),
  xpos(0*cm),
  ypos(0*cm),
  zpos(0*cm),
  _sciTilt(0),
  _sciWidth(0.1*cm),
  _tungstenWidth(0.1*cm),
  _unmag(1*cm),
  _sciNum(-1),
  _region(NULL),
  active(0),
  absorberactive(0),
  layer(lyr),
  undulations(5),
  usedrawing(0)
{}

//_______________________________________________________________
//_______________________________________________________________
int PHG4AccDetector::IsInCylinderActive(const G4VPhysicalVolume * volume)
{
  //  cout << "checking detector" << endl;
  if (active && scinti_vol.find(volume) != scinti_vol.end())
    {
      return scinti_vol.find(volume)->second;
    }
  if (absorberactive && volume == cylinder_physi)
    {
      return -1;
    }
  return INACTIVE;
}


//_______________________________________________________________
void PHG4AccDetector::Construct( G4LogicalVolume* logicWorld )
{
  TrackerMaterial = G4Material::GetMaterial(material);

  if ( ! TrackerMaterial)
    {
      std::cout << "Error: Can not set material" << std::endl;
      exit(-1);
    }
  if (usedrawing)
    {
      // hardcoded - 1mm scintillator, 2 mm tungsten
      cout << "using geometry from engineering drawing:" << endl;
      cout << "Scintillator width 1mm" << endl;
      cout << "Tungsten width 2mm" << endl;
      cout << "Thickness 10cm" << endl;
      _sciWidth = 1*mm;
      _tungstenWidth = 2*mm;
      TrackerThickness = 10*cm; // this is already converted to internal G4 units
      G4VisAttributes* scintillatorVis = new G4VisAttributes();
      G4VisAttributes* CylinderVis = new G4VisAttributes();
      scintillatorVis->SetVisibility(true);
      scintillatorVis->SetForceSolid(true);
      PHG4Utils::SetColour(scintillatorVis, "G4_POLYSTYRENE");
      CylinderVis->SetVisibility(true);
      CylinderVis->SetForceSolid(false);
      PHG4Utils::SetColour(CylinderVis, material);
      G4Tubs* _cylinder_solid = new G4Tubs(G4String(GetName().c_str()),
					   radius,
					   radius + TrackerThickness,
					   length / 2.0, 0, twopi);
      cylinder_logic = new G4LogicalVolume(_cylinder_solid,
					   TrackerMaterial,
					   G4String(GetName().c_str()),
					   0, 0, 0);
      cylinder_logic->SetVisAttributes(CylinderVis);
      cylinder_physi = new G4PVPlacement(0, G4ThreeVector(xpos, ypos, zpos),
					 cylinder_logic,
					 G4String(GetName().c_str()),
					 logicWorld, 0, false, overlapcheck);


      G4Box *box1 = new  G4Box("TMPBOX", 100.0*mm / 2., 6.66*mm / 2., (length / 2.0)*1.1);
      double phimin = 3 * M_PI / 2. - 19.47 / 180.*M_PI;
      double phimin2 = M_PI / 2. - 19.47 / 180.*M_PI;
      double phiall = 2 * 19.47 / 180.*M_PI;
      double phiallhalf = 19.47 / 180.*M_PI;

      // this one is _/
      G4Tubs *cylcut0 =   new G4Tubs("LOWERHALFARC",
				     5*cm,
				     5*cm + 1*mm,
				     length / 2.0, 3*M_PI / 2. , phiallhalf);
      /* this one is /-\  (and /-\ inside a // comment makes the compiler barf) */
	   G4Tubs *cylcut1 =   new G4Tubs("UPPERFULLARC",
					  5*cm,
					  5*cm + 1*mm,
					  length / 2.0, phimin2, phiall);
      // this one is \_/
      G4Tubs *cylcut2 =   new G4Tubs("LOWERFULLARC",
				     5*cm,
				     5*cm + 1*mm,
				     length / 2.0, phimin, phiall);
      // this one is /-
      G4Tubs *cylcut3 =   new G4Tubs("UPPERHALFARC",
				     5*cm,
				     5*cm + 1*mm,
				     length / 2.0, M_PI / 2. , phiallhalf);
      // this one is just to move the scintillator arc to zero, if we don't do this
      // the arc will sit at 5cm with respect to zero and the later positioning
      // of the combined volume will be a lot more complicated
      double tenbythree = 100./3.;
      G4ThreeVector xshift0(0, 5*cm + 1*mm - 6.66 / 2.*mm, 0);
      G4IntersectionSolid *movezero = new G4IntersectionSolid("ZEROIT", box1, cylcut0, 0, xshift0);
      G4ThreeVector xshift1(tenbythree*mm - 0.5*mm, -5*cm - 1*mm + 6.66 / 2.*mm, 0);
      G4UnionSolid *zeroone =  new G4UnionSolid("ZEROONE", movezero, cylcut1, 0, xshift1);
      G4ThreeVector xshift2(2*(tenbythree*mm - 0.5*mm), 5*cm + 1*mm - 6.66 / 2.*mm, 0);
      G4UnionSolid *zeroonetwo = new G4UnionSolid("ZEROONETWO", zeroone, cylcut2, 0, xshift2);
      G4ThreeVector xshift3(3*(tenbythree*mm - 0.5*mm), -5*cm - 1*mm + 6.66 / 2.*mm, 0);
      G4UnionSolid *zeroonetwothree = new G4UnionSolid("ZEROONETWOTHREE", zeroonetwo, cylcut3, 0, xshift3);
      G4LogicalVolume *cylinder_logic4 = new G4LogicalVolume(zeroonetwothree, G4Material::GetMaterial("G4_POLYSTYRENE"), "SCINTILLATOR", 0, 0, 0);
      cylinder_logic4->SetVisAttributes(scintillatorVis);
      if (_sciNum <= 0)
        {
          double circumference = radius * 2 * M_PI;

          int    sciNumsci = int(circumference / (_sciWidth + _tungstenWidth));
          _sciNum = sciNumsci;
        }
      cout << "using " << _sciNum << " scintillator slats" << endl;
      double phi_increment = 360. / _sciNum;
      ostringstream slatname;
      double Rmid = radius;
      for (int i = 0; i < _sciNum; i++)
        {
          double phi = i  * phi_increment;
          G4ThreeVector myTrans = G4ThreeVector(Rmid * cos(phi * M_PI / 180.), Rmid * sin(phi * M_PI / 180.), 0);
          G4RotationMatrix Rot(0, 0, 0);
          Rot.rotateZ(phi*deg);
          slatname.str("");
          slatname << "SLAT_" << i;
          G4VPhysicalVolume* box_vol_tmp = new G4PVPlacement(G4Transform3D(Rot, G4ThreeVector(myTrans)),
							     cylinder_logic4,
							     G4String(slatname.str()),
							     cylinder_logic, 0, false, overlapcheck);
          scinti_vol[box_vol_tmp] = i;
        }
    }
  else
    {
  G4VisAttributes* scintillatorVis = new G4VisAttributes();
  G4VisAttributes* CylinderVis = new G4VisAttributes();
  scintillatorVis->SetVisibility(true);
  scintillatorVis->SetForceSolid(true);
  PHG4Utils::SetColour(scintillatorVis, "G4_POLYSTYRENE");
  CylinderVis->SetVisibility(true);
  CylinderVis->SetForceSolid(true);
  PHG4Utils::SetColour(CylinderVis, material);
  G4Tubs* _cylinder_solid = new G4Tubs(G4String(GetName().c_str()),
                                       radius,
                                       radius + TrackerThickness,
                                       length / 2.0, 0, twopi);
  cylinder_logic = new G4LogicalVolume(_cylinder_solid,
                                       TrackerMaterial,
                                       G4String(GetName().c_str()),
                                       0, 0, 0);
  cylinder_logic->SetVisAttributes(CylinderVis);
  cylinder_physi = new G4PVPlacement(0, G4ThreeVector(xpos, ypos, zpos),
                                     cylinder_logic,
                                     G4String(GetName().c_str()),
                                     logicWorld, 0, false, overlapcheck);


  // now for the undulating slats:
  // I could not find a way to just use a sin/cos to create this volume. So here
  // I start with a solid elliptical tube, create a smaller tube and subtract it so I am left
  // with a hollow elliptical tube. Then I create a box and subtract that one which leaves me
  // with a left and right half elliptical tube. Then I add a few of those vertically which
  // results in the desired undulating slat

  // if number of scintillators not set, divide up our cylinder
  if (_sciNum <= 0)
    {
      double circumference = radius * 2 * M_PI;

      int    sciNumsci = int(circumference / (_sciWidth + _tungstenWidth));
      _sciNum = sciNumsci;
    }
  double ell_y = ((TrackerThickness + (undulations - 1) * _sciWidth) / undulations) / 2.;
  double ell_x = _unmag / 2.0;
  double ell_len = length / 2.0;
  ell = new G4EllipticalTube("ELL", ell_x, ell_y, ell_len);
  G4EllipticalTube* ell_inner = new G4EllipticalTube("ELL", ell_x - _sciWidth, ell_y - _sciWidth, ell_len*1.1);
  G4ThreeVector zero(0, 0, 0);
  G4SubtractionSolid* subtract_ell = new G4SubtractionSolid("ELL_HOLLOW", ell, ell_inner, 0, zero);
  G4Box *box = new G4Box("BOX", ell_x, ell_y, ell_len*1.1);
  G4ThreeVector leftshift(-ell_x, 0, 0);
  G4ThreeVector rightshift(ell_x, 0, 0);
  G4SubtractionSolid* right_ell = new G4SubtractionSolid("RIGHT_ELL", subtract_ell, box, 0, leftshift);
  G4SubtractionSolid* left_ell = new G4SubtractionSolid("LEFT_ELL", subtract_ell, box, 0, rightshift);
  G4ThreeVector upshift_0(0, 2*ell_y - _sciWidth, 0);
  G4ThreeVector upshift(0, 2.*(2*ell_y - _sciWidth), 0);
  G4ThreeVector upshift_2(0, 4.*(2*ell_y - _sciWidth), 0);
  G4UnionSolid *lr_ell_0 = new G4UnionSolid("LR_ELL_0", left_ell, right_ell, 0, upshift_0);
  G4UnionSolid *lr_ell_2 = new G4UnionSolid("LR_ELL_2", lr_ell_0, lr_ell_0, 0, upshift);
  G4UnionSolid *lr_ell = new G4UnionSolid("LR_ELL", lr_ell_2, left_ell, 0, upshift_2);
  G4Material* accordion_mat = G4Material::GetMaterial("G4_POLYSTYRENE");
  scinti_logic = new G4LogicalVolume(lr_ell,
                                     accordion_mat,
                                     G4String(GetName().c_str()),
                                     0, 0, 0);
  scinti_logic->SetVisAttributes(scintillatorVis);
  //   G4VPhysicalVolume* ell_vol_tmp = new G4PVPlacement(0, G4ThreeVector(xpos, ypos, zpos),
  // 						     cylinder_logic,
  // 						     G4String(GetName().c_str()),
  // 						     logicWorld, 0, false, overlapcheck);
  double phi_increment = 360. / _sciNum;
  ostringstream slatname;
  double Rmid = radius + ell_y;
  for (int i = 0; i < _sciNum; i++)
    {
      double phi = i  * phi_increment;
      G4ThreeVector myTrans = G4ThreeVector(Rmid * -1.*sin(phi * M_PI / 180.), Rmid * cos(phi * M_PI / 180.), 0);
      G4RotationMatrix Rot(0, 0, 0);
      Rot.rotateZ(phi*deg);
      slatname.str("");
      slatname << "SLAT_" << i;
      G4VPhysicalVolume* box_vol_tmp = new G4PVPlacement(G4Transform3D(Rot, G4ThreeVector(myTrans)),
							 scinti_logic,
							 G4String(slatname.str()),
							 cylinder_logic, 0, false, overlapcheck);
      scinti_vol[box_vol_tmp] = i;
    }

    }

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
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeomv2(radius / cm, (zpos - length / 2.) / cm, (zpos + length / 2.) / cm, TrackerThickness / cm, _sciNum );
      geo->AddLayerGeom(layer, mygeom);
    }
  return;
}

void
PHG4AccDetector::Print(const std::string &what) const
{
  cout << "radius: " << radius << endl;
  return;
}

