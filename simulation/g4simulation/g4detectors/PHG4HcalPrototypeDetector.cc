// This is G4 detector implementation for the hcal prototype
// Created on 1/27/2014, Liang, HeXC
// Updated on 3/21/2014, Liang, HeXC
// Updated on 4/18/2014, Liang, HeXC, Updating detector construction (proper tilting angles!)

#include "PHG4HcalPrototypeDetector.h"

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Material.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Trap.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <Geant4/G4NistManager.hh>
#include <Geant4/G4GeometryManager.hh>
#include <Geant4/G4PhysicalVolumeStore.hh>
#include <Geant4/G4LogicalVolumeStore.hh>
#include <Geant4/G4SolidStore.hh>

#include <Geant4/G4RunManager.hh>

#include <sstream>

using namespace std;

// static double no_overlap = 0.00015 * cm; // added safety margin against overlaps by using same boundary between volumes
static const double inch = 2.54 * cm;

PHG4HcalPrototypeDetector::PHG4HcalPrototypeDetector( PHCompositeNode *Node, const std::string &dnam, const int lyr  ):
  PHG4Detector(Node, dnam),
  physiWorld(nullptr),
  logicWorld(nullptr),
  logicHcalBox(nullptr),
  solidHcalBox(nullptr),
  physiHcalBox(nullptr),
  logicHcal2ScintLayer(nullptr),
  logicHcal1ScintLayer(nullptr),
  logicHcal2AbsLayer(nullptr),
  logicHcal1AbsLayer(nullptr),
  world_mat(nullptr),
  steel(nullptr),
  scint_mat(nullptr),
  outer_tilt_angle(NAN),
  active(0),
  absorberactive(0),
  layer(lyr),
  blackhole(0)
{
  // Outer Hcal dimensions
  hcal2ScintSizeX = 2.54*26.1*cm;      // from Don's drawing
  hcal2ScintSizeZ = 2.54*29.53*cm;     // from Don's drawing
  hcal2ScintSizeY = 0.85*cm;           // Changed from 8.0cm to 8.5cm following Don's suggestion

  // Inner Hcal dimensions
  hcal1ScintSizeX = 2.54*12.43*cm;     // from Don's drawing
  hcal1ScintSizeZ = 2.54*17.1*cm;      // from Don's drawing
  hcal1ScintSizeY = 0.85*cm;           // Changed from 8.0cm to 8.5cm following Don's suggestion

  // Radius dimension
  hcal1RadiusIn = 1855*mm;
  hcal2RadiusIn = hcal1RadiusIn + hcal1ScintSizeX;

  nScint360 = 15;  // This is legacy variable in case one wants to build a cylindrical calorimeter 
  nHcal2Layers = nHcal1Layers = 16;

  hcal2DPhi = 0.288;  // calculated based on Don's drawing.  (twopi/16.0);
  hcal1DPhi = 0.27;   // calculated based on Don's drawing.  (twopi/16.0);

  hcal2TiltAngle = 0.14;
  hcal1TiltAngle = 0.27;

  hcal2Abs_dxa = 2.54*27.1*cm;   // these numbers are from Don't drawings
  hcal2Abs_dxb = hcal2Abs_dxa;
  hcal2Abs_dya = 2.54*1.099*cm;
  hcal2Abs_dyb = 2.54*1.776*cm;  
  hcal2Abs_dz  = hcal2ScintSizeZ;

  hcalJunctionSizeX = 2.54*1.0*cm;
  hcalJunctionSizeY = hcal2ScintSizeY;
  hcalJunctionSizeZ = hcal2Abs_dz;

  hcal1Abs_dxa = hcal1ScintSizeX;  // these numbers are from Don't drawings
  hcal1Abs_dxb = hcal1ScintSizeX;
  hcal1Abs_dya = 2.54*0.787*cm;
  hcal1Abs_dyb = 2.54*1.115*cm;    
  hcal1Abs_dz  = hcal1ScintSizeZ;

  // The frame box for the HCAL.
  hcalBoxSizeX = 1.1*(hcal2ScintSizeX + hcal1ScintSizeX); 
  hcalBoxSizeY = 2.54*40.0*cm;    // need to get more accurate dimension for the box             
  hcalBoxSizeZ = 1.1*hcal2ScintSizeZ;

  hcalBoxRotationAngle_z = 0.0*rad;   // Rotation along z-axis
  hcalBoxRotationAngle_y = 0.0*rad;   // Rotation along y-axis

  // create commands for interactive definition of the detector  
  fDetectorMessenger = new PHG4HcalPrototypeDetectorMessenger(this);
}

//_______________________________________________________________
//_______________________________________________________________
int
PHG4HcalPrototypeDetector::IsInHcalPrototype(G4VPhysicalVolume * volume) const
{
  /*
  if (active && volume == sandwich_vol[2])
    {
      return 1;
    }
  if (absorberactive && (volume == sandwich_vol[0] || volume == sandwich_vol[1]))
    {
      return -1;
    }
  if (absorberactive && sandwich_vol[3] && volume == sandwich_vol[3])
    {
      return -1;
    }
  */

  return 0;  // not sure what value to return for now.
}

void PHG4HcalPrototypeDetector::Construct( G4LogicalVolume* world )
{
  logicWorld = world;

  //  physiWorld = world->GetSolid();

  DefineMaterials();
  
  ConstructDetector();

  /*

  hcal_rotm.rotateX(x_rot);
  hcal_rotm.rotateY(y_rot);
  hcal_rotm.rotateZ(z_rot);

  new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(place_in_x, place_in_y, place_in_z)), hcal_log, "Hcal", logicWorld, 0, false, overlapcheck);
  ConstructOuterHcal(hcal_log);
    return;
  */
}


void PHG4HcalPrototypeDetector::DefineMaterials()
{ 

  // Water is defined from NIST material database
  G4NistManager * nist = G4NistManager::Instance();
  
  world_mat = nist->FindOrBuildMaterial("G4_AIR");

  steel = nist->FindOrBuildMaterial("G4_Fe"); // Need to fix this *******
  scint_mat = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  
  // If one wishes to test other density value for water material, one should use instead:
  //G4Material * H2O = man->BuildMaterialWithNewDensity("G4_WATER_MODIFIED","G4_WATER",1000*g/cm/cm/cm);
  //G4cout << "-> Density of water material (g/cm3)=" << waterMaterial->GetDensity()/(g/cm/cm/cm) << G4endl;
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

// We now build our detector solids, etc
// 
G4VPhysicalVolume*  PHG4HcalPrototypeDetector::ConstructDetector()
{

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // HCAL Frame Box for enclosing the inner and the outer hcal
  // which allows us to rotate the whole calorimeter easily
  solidHcalBox = new G4Box("HcalBox",				           //its name
			 hcalBoxSizeX/2,hcalBoxSizeY/2,hcalBoxSizeZ/2);    //its size

  logicHcalBox = new G4LogicalVolume(solidHcalBox,      //its solid
				     world_mat,	        //its material
				     "HcalBox");        //its name

  // Place the HcalBox inside the World
  G4ThreeVector boxVector = G4ThreeVector(hcal1RadiusIn + hcalBoxSizeX/2.0,  0, 0); 
  G4RotationMatrix rotBox  = G4RotationMatrix();

  rotBox.rotateZ(hcalBoxRotationAngle_z);     // Rotate the whole calorimeter box along z-axis
  rotBox.rotateY(hcalBoxRotationAngle_y);     // Rotate the whole calorimeter box along z-axis

  G4Transform3D boxTransform = G4Transform3D(rotBox, boxVector);
  
  physiHcalBox = new G4PVPlacement(boxTransform,        // rotation + positioning
				   logicHcalBox,	//its logical volume
				   "HcalBox",		//its name
				   logicWorld,		//its mother  volume
				   false,		//no boolean operation
				   0,			//copy number
				   checkOverlaps);      //checking overlaps
  
  // HCal Box visualization attributes
  G4VisAttributes* hcalBoxVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0)); // blue
  hcalBoxVisAtt->SetVisibility(true);
  //hcalBoxVisAtt->SetForceSolid(true);
  logicHcalBox->SetVisAttributes(hcalBoxVisAtt);


  // Make outer hcal section
  // 
  // Build scintillator layer box as a place holder for the scintillator sheets
  // Make it 0.05cm thinner than the true gap width for avoiding geometry overlaps
  G4Box* hcal2ScintLayer = new G4Box("hcal2ScintLayer",		                              //its name
				     hcal2ScintSizeX/2,(hcal2ScintSizeY-0.05*cm)/2,hcal2ScintSizeZ/2);  //its size
  
  
  logicHcal2ScintLayer = 
    new G4LogicalVolume(hcal2ScintLayer,      // its solid name
			world_mat,            // material, the same as the material of the world valume
			"hcal2ScintLayer");   // its name
  
  // Scintillator visualization attributes
  G4VisAttributes* scintVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //Red
  scintVisAtt->SetVisibility(true);
  //scintVisAtt->SetForceSolid(true);
  logicHcal2ScintLayer->SetVisAttributes(scintVisAtt);

  // Construct outer scintillator 1U
  G4double outer1UpDz = 649.8*mm;
  //  G4double outer1UpTheta = 0.0;
  //G4double outer1UpPhi = 0.0;
  G4double outer1UpDy1 = 2*0.35*cm;
  //G4double outer1UpDy2 = 0.35*cm;
  //G4double outer1UpDx1 = 0.25*179.3*mm;
  G4double outer1UpDx2 = 179.3*mm;
  //G4double outer1UpDx3 = 0.25*113.2*mm;
  G4double outer1UpDx4 = 113.2*mm;
  //G4double outer1UpAlp1 = 0., outer1UpAlp2 = 0.;

  G4Trap *outer1USheetSolid = new G4Trap("outer1USheet",                      //its name
					 outer1UpDy1, outer1UpDz, 
					 outer1UpDx2, outer1UpDx4); //its size
    
  G4LogicalVolume *logicOuter1USheet = new G4LogicalVolume(outer1USheetSolid, 
							 scint_mat,
							 "outer1USheet");

  // Scintillator visualization attributes
  G4VisAttributes* scintSheetVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.0)); //White
  scintSheetVisAtt->SetVisibility(true);
  scintSheetVisAtt->SetForceSolid(true);

  logicOuter1USheet->SetVisAttributes(scintSheetVisAtt);

  // Detector placement
  // Position the outer right 1U sheet
  G4ThreeVector threeVecOuter1U_1 = G4ThreeVector(0*cm, 0*cm, -0.252*(outer1UpDx2+outer1UpDx4));
  G4RotationMatrix rot1U_1  = G4RotationMatrix();
  rot1U_1.rotateZ(90*deg);
  rot1U_1.rotateX(-90*deg);

  G4Transform3D transformOuter1U_1 = G4Transform3D(rot1U_1,threeVecOuter1U_1);

  new G4PVPlacement(transformOuter1U_1,
		    logicOuter1USheet,
		    "outer1USheet",
		    logicHcal2ScintLayer,
		    false,
		    0,                    // Copy one
		    checkOverlaps);

  
  // Detector placement
  // Position the outer left 1U sheet
  G4ThreeVector threeVecOuter1U_2 = G4ThreeVector(0*cm, 0*cm, 0.252*(outer1UpDx2+outer1UpDx4));
  G4RotationMatrix rot1U_2  = G4RotationMatrix();
  rot1U_2.rotateZ(90*deg);
  rot1U_2.rotateX(90*deg);

  G4Transform3D transformOuter1U_2 = G4Transform3D(rot1U_2,threeVecOuter1U_2);

  new G4PVPlacement(transformOuter1U_2,
		    logicOuter1USheet,
		    "outer1USheet",
		    logicHcal2ScintLayer,
		    false,
		    1,                      // Copy two
		    checkOverlaps);

  // Construct outer scintillator 2U 
  G4double outer2UpDz = 0.5*649.8*mm;
  G4double outer2UpTheta = 8.8*M_PI/180.;
  G4double outer2UpPhi = 0.0*M_PI/180.;
  G4double outer2UpDy1 = 0.35*cm;
  G4double outer2UpDy2 = 0.35*cm;
  G4double outer2UpDx1 = 0.5*179.3*mm, outer2UpDx2 = 0.5*179.3*mm;
  G4double outer2UpDx3 = 0.5*113.2*mm, outer2UpDx4 = 0.5*113.2*mm;
  G4double outer2UpAlp1 = 0.*M_PI/180., outer2UpAlp2 = 0.*M_PI/180;

  G4Trap *outer2USheetSolid = new G4Trap("outer2USheet",
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

  G4LogicalVolume *logicOuter2USheet = new G4LogicalVolume(outer2USheetSolid, 
							 scint_mat,
							 "outer2USheet");

  logicOuter2USheet->SetVisAttributes(scintSheetVisAtt);

  
  // Detector placement
  // Position the right most sheet
  G4ThreeVector threeVecOuter2U_1 = G4ThreeVector(0*cm, 0*cm, -0.755*(outer1UpDx2+outer1UpDx4));
  G4RotationMatrix rot2U_1  = G4RotationMatrix();
  rot2U_1.rotateY(-90*deg);

  G4Transform3D transformOuter2U_1 = G4Transform3D(rot2U_1,threeVecOuter2U_1);

  new G4PVPlacement(transformOuter2U_1,
		    logicOuter2USheet,
		    "outer2USheet",
		    logicHcal2ScintLayer,
		    false,
		    0,                     // copy one
		    checkOverlaps);

  
  // Detector placement
  // Position the outer left most sheet
  G4ThreeVector threeVecOuter2U_2 = G4ThreeVector(0*cm, 0*cm, 0.755*(outer1UpDx2+outer1UpDx4));
  G4RotationMatrix rot2U_2  = G4RotationMatrix();
  rot2U_2.rotateY(-90*deg);
  rot2U_2.rotateX(180*deg);

  G4Transform3D transformOuter2U_2 = G4Transform3D(rot2U_2,threeVecOuter2U_2);

  new G4PVPlacement(transformOuter2U_2,
		    logicOuter2USheet,
		    "outer2USheet",
		    logicHcal2ScintLayer,
		    false,
		    1,                     // copy two
		    checkOverlaps);  
  
  G4Trap* hcal2AbsLayer =    
    new G4Trap("hcal2AbsLayer",                      //its name
	       hcal2Abs_dz, hcal2Abs_dxa, 
	       hcal2Abs_dyb, hcal2Abs_dya); //its size
  
  // Add extra absorber materials at the junction
  G4Box* hcalJunction = new G4Box("HcalJunction",		     //its name
				  hcalJunctionSizeX/2, hcalJunctionSizeY/2, hcalJunctionSizeZ/2);  //its size
  
  // define the transformation for placing the junction to the inner side of the hcal2 absorber layer
  G4ThreeVector threeVecJunction = G4ThreeVector( -hcal2Abs_dya/2.0 - 1.1*hcalJunctionSizeY , hcal2Abs_dxa/2.0 - hcalJunctionSizeX/2, 0.0*cm);
  G4RotationMatrix rotJunction  = G4RotationMatrix();
  rotJunction.rotateZ(90*deg);
  rotJunction.rotateX(0*deg);
  G4Transform3D transformJunction = G4Transform3D(rotJunction,threeVecJunction);
  
  // attach the junction to the absorber using G4UnionSolid
  G4UnionSolid* hcal2AbsLayerJunct = new G4UnionSolid("hcal2AbsLayerJunct", hcal2AbsLayer, hcalJunction, transformJunction);
  
  logicHcal2AbsLayer = new G4LogicalVolume(hcal2AbsLayerJunct,    //  hcal2AbsLayer,
					   steel,
					   "hcal2AbsLayer");
  /* old plain absorber
  logicHcal2AbsLayer = new G4LogicalVolume(hcal2AbsLayer,
					   steel,
					   "hcal2AbsLayer");
  */

  // hcal2Aborber visualization attributes
  G4VisAttributes* hcal2AbsVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,1.0)); 
  hcal2AbsVisAtt->SetVisibility(true);
  logicHcal2AbsLayer->SetVisAttributes(hcal2AbsVisAtt);
      
  // place scintillator layers insides the HCal2 absorber volume 
  //
  
  // Calculate the title angle
  // SetTiltViaNcross(2);    // with 2 crossings?  This function does not work yet!  :((

  G4double theta = 0.;
  G4double theta2 = hcal2TiltAngle;   // I am not convinced that I got the tilt angle right.  So I simply forced it to 
  G4double rScintLayerCenter = hcal2RadiusIn + hcal2ScintSizeX/2.0 + 2.54*1.0*cm;  // move the scintillator toward the rear end of HCAL2
  G4double RmidScintLayerX =  0.81*(hcal2ScintSizeX - hcal1ScintSizeX)/2.0 + 2.54*1.0*cm + 2*2.54*cm;  // move the scintillator toward the rear end of HCAL2
  G4double rAbsLayerCenter = hcal2RadiusIn + hcal2Abs_dxa/2.0;
  G4double RmidAbsLayerX =  0.81*(hcal2Abs_dxa - hcal1ScintSizeX)/2.0 + 2*2.54*cm;
  G4int LayerNum = 1;
  G4double xPadding = 0.0;
  G4double yPadding = 0.0;
  G4double tiltPadding = 0.0;

  for (G4int iLayer = -nHcal2Layers/2; iLayer < nHcal2Layers/2; iLayer++ ) { 
    //
    // Place outer absorber layer first
    //
    theta = hcal2DPhi/nHcal2Layers * iLayer;
    G4double xposShift = rAbsLayerCenter * (cos(theta) - 1.0);
    G4double yposShift = rAbsLayerCenter * sin(theta);
    xPadding = 0.013*RmidAbsLayerX*LayerNum;      // adding an extra padding for placing the absorber and scintillator 
    G4ThreeVector absTrans = G4ThreeVector(RmidAbsLayerX + xposShift - xPadding, yposShift, 0); 
    G4RotationMatrix rotAbsLayer  = G4RotationMatrix();
    rotAbsLayer.rotateY(90*deg);
    rotAbsLayer.rotateX(90*deg);
    rotAbsLayer.rotateY(-90*deg);
    tiltPadding = LayerNum*0.005*rad;
    rotAbsLayer.rotateZ((iLayer + nHcal2Layers/2)*0.02*rad + tiltPadding);   // Added a funge increment factor 0.01
    
    G4Transform3D transformAbs = G4Transform3D(rotAbsLayer,absTrans);
    
    new G4PVPlacement(transformAbs,            //rotation,position
		      logicHcal2AbsLayer,      //its logical volume
		      "hcal2AbsLayer",         //its name
		      logicHcalBox,            //its mother  volume
		      false,                   //no boolean operation
		      LayerNum - 1,            //copy number
		      checkOverlaps);          //overlaps checking 
    //
    // Place outer hcal scintillator layer
    //
    theta += 0.5*hcal2DPhi/nHcal2Layers; 
    G4cout << "M_PI: " << M_PI << "     theta: " << theta << "    TileAngle: " << theta2 << G4endl;
    xposShift = rScintLayerCenter * (cos(theta) - 1.0);
    yposShift = rScintLayerCenter * sin(theta);
    G4ThreeVector myTrans = G4ThreeVector(RmidScintLayerX + xposShift - xPadding, yposShift+0.3*cm, 0); // hcal2RadiusIn + hcal1ScintSizeX/2.0
    G4RotationMatrix rotm  = G4RotationMatrix();
    //    rotm.rotateZ((theta+theta2+0.02)*rad);        // Added a funge increment factor 0.02
    rotm.rotateZ((iLayer + nHcal2Layers/2 + 1)*0.020*rad+tiltPadding);        // Added a funge increment factor 0.02
    G4Transform3D transform = G4Transform3D(rotm,myTrans);
     
    G4cout << "  iLayer " << iLayer << G4endl;

    new G4PVPlacement(transform,             //rotation,position
                      logicHcal2ScintLayer,  //its logical volume
                      "hcal2ScintLayer",     //its name
                      logicHcalBox,          //its mother volume
                      false,                 //no boolean operation
                      LayerNum - 1,          //copy number
                      checkOverlaps);        //checking overlaps 
    LayerNum++;

  }

  // Make INNER hcal section
  // 
  // Build scintillator layer box as a place holder for the scintillator sheets
  // Make it 0.05cm thinner than the true gap width for avoiding geometry overlaps
  G4Box* hcal1ScintLayer = new G4Box("hcal1ScintLayer",		                              //its name
				     hcal1ScintSizeX/2,(hcal1ScintSizeY-0.05*cm)/2,hcal1ScintSizeZ/2);  //its size

  logicHcal1ScintLayer = 
    new G4LogicalVolume(hcal1ScintLayer,     // its solid name
			world_mat,           // simply use world material
			"hcal1ScintLayer");  // its name

  logicHcal1ScintLayer->SetVisAttributes(scintVisAtt);

  // Constructing scintillator sheets
  // Construct inner scintillator 1U
  // The numbers are read off from Don's drawings
  G4double inner1UpDz = 316.8*mm;
  //  G4double inner1UpTheta = 0.0;
  //G4double inner1UpPhi = 0.0;
  G4double inner1UpDy1 = 2*0.35*cm;
  //G4double inner1UpDy2 = 0.35*cm;
  //G4double inner1UpDx1 = 0.25*108.6*mm;
  G4double inner1UpDx2 = 108.6*mm;
  //G4double inner1UpDx3 = 0.25*77.4*mm;
  G4double inner1UpDx4 = 77.4*mm;
  //G4double inner1UpAlp1 = 0., inner1UpAlp2 = 0.;

  G4Trap *inner1USheetSolid = new G4Trap("inner1USheet",             //its name
					 inner1UpDy1, inner1UpDz, 
					 inner1UpDx2, inner1UpDx4);  //its size
    
  G4LogicalVolume *logicInner1USheet = new G4LogicalVolume(inner1USheetSolid, 
							 scint_mat,
							 "inner1USheet");
  logicInner1USheet->SetVisAttributes(scintSheetVisAtt);

  // Detector placement
  // Position the inner right 1U sheet
  G4ThreeVector threeVecInner1U_1_inner = G4ThreeVector(0*cm, 0*cm, -0.252*(inner1UpDx2+inner1UpDx4));
  //G4RotationMatrix rot1U_1  = G4RotationMatrix();   // We could reuse this from outer hcal construct
  //rot1U_1.rotateZ(90*deg);
  //rot1U_1.rotateX(-90*deg);

  G4Transform3D transformInner1U_1 = G4Transform3D(rot1U_1,threeVecInner1U_1_inner);

  new G4PVPlacement(transformInner1U_1,
		    logicInner1USheet,
		    "inner1USheet",
		    //logicWorld,
		    logicHcal1ScintLayer,
		    false,
		    0,                    // Copy one
		    checkOverlaps);

  
  // Detector placement
  // Position the inner left 1U sheet
  G4ThreeVector threeVecInner1U_2_inner = G4ThreeVector(0*cm, 0*cm, 0.252*(inner1UpDx2+inner1UpDx4));
  //G4RotationMatrix rot1U_2  = G4RotationMatrix();
  //rot1U_2.rotateZ(90*deg);
  //rot1U_2.rotateX(90*deg);

  G4Transform3D transformInner1U_2 = G4Transform3D(rot1U_2,threeVecInner1U_2_inner);

  new G4PVPlacement(transformInner1U_2,
		    logicInner1USheet,
		    "inner1USheet",
		    //logicWorld,
		    logicHcal1ScintLayer,
		    false,
		    1,                      // Copy two
		    checkOverlaps);

  // Construct inner scintillator 2U 
  G4double inner2UpDz = 0.5*316.8*mm;
  G4double inner2UpTheta = 8.8*M_PI/180.;
  G4double inner2UpPhi = 0.0*M_PI/180.;
  G4double inner2UpDy1 = 0.35*cm;
  G4double inner2UpDy2 = 0.35*cm;
  G4double inner2UpDx1 = 0.5*108.6*mm, inner2UpDx2 = 0.5*108.6*mm;
  G4double inner2UpDx3 = 0.5*77.4*mm, inner2UpDx4 = 0.5*77.4*mm;
  G4double inner2UpAlp1 = 0.*M_PI/180., inner2UpAlp2 = 0.*M_PI/180;

  G4Trap *inner2USheetSolid = new G4Trap("inner2USheet",
				  inner2UpDz,
				  inner2UpTheta,
				  inner2UpPhi, 
				  inner2UpDy1,
				  inner2UpDx1,
				  inner2UpDx2,
				  inner2UpAlp1,
				  inner2UpDy2,
				  inner2UpDx3,
				  inner2UpDx4,
				  inner2UpAlp2);

  G4LogicalVolume *logicInner2USheet = new G4LogicalVolume(inner2USheetSolid, 
							 scint_mat,
							 "inner2USheet");

  logicInner2USheet->SetVisAttributes(scintSheetVisAtt);

  
  // Detector placement
  // Position the right most sheet
  G4ThreeVector threeVecInner2U_1_inner = G4ThreeVector(0*cm, 0*cm, -0.755*(inner1UpDx2+inner1UpDx4));
  //G4RotationMatrix rot2U_1  = G4RotationMatrix();
  //rot2U_1.rotateY(-90*deg);

  G4Transform3D transformInner2U_1 = G4Transform3D(rot2U_1,threeVecInner2U_1_inner);

  new G4PVPlacement(transformInner2U_1,
		    logicInner2USheet,
		    "inner2USheet",
		    //logicWorld,
		    logicHcal1ScintLayer,
		    false,
		    0,                     // copy one
		    checkOverlaps);

  
  // Detector placement
  // Position the inner left most sheet
  G4ThreeVector threeVecInner2U_2_inner = G4ThreeVector(0*cm, 0*cm, 0.755*(inner1UpDx2+inner1UpDx4));
  //G4RotationMatrix rot2U_2  = G4RotationMatrix();
  //rot2U_2.rotateY(-90*deg);
  //rot2U_2.rotateX(180*deg);

  G4Transform3D transformInner2U_2 = G4Transform3D(rot2U_2,threeVecInner2U_2_inner);

  new G4PVPlacement(transformInner2U_2,
		    logicInner2USheet,
		    "inner2USheet",
		    //logicWorld,
		    logicHcal1ScintLayer,
		    false,
		    1,                     // copy two
		    checkOverlaps);  
  

  // Build hcal1 absorber layer in trapezoid shape 
  //    
  /*
  G4double hcal1Abs_dxa = hcal1ScintSizeX; // hcal1Abs_dxb = hcal1ScintSizeX;
  G4double hcal1Abs_dya = 0.787*2.54*cm, hcal1Abs_dyb = 1.115*2.54*cm;    
  // these numbers are from Don't drawings
  G4double hcal1Abs_dz  = hcal1ScintSizeZ; 
  */
  G4Trap* hcal1AbsLayer =    
    new G4Trap("hcal1AbsLayer",                 //its name
	       hcal1Abs_dz, hcal1Abs_dxa, 
	       hcal1Abs_dyb, hcal1Abs_dya);      //its size
  
  logicHcal1AbsLayer = new G4LogicalVolume(hcal1AbsLayer,
					   steel,
							 "hcal1AbsLayer"); 

  logicHcal1AbsLayer->SetVisAttributes(hcal2AbsVisAtt);

  // place scintillator layers insides the HCal1 absorber volume 
  //
  // Calculate the title angle

  theta2 = hcal1TiltAngle;
  rAbsLayerCenter = hcal1RadiusIn + hcal1ScintSizeX / 2.0;
  //  RmidAbsLayerX =  (hcal2Abs_dxa - hcal1ScintSizeX)/2.0;
  RmidAbsLayerX =  hcal2Abs_dxa/2.0;
  LayerNum = 1;          // These three parameters are for making fine adjustment of the placement of inner hcal
  yPadding = -2.54*0.8*cm;
  for (G4int iLayer = -nHcal2Layers/2 -1; iLayer < nHcal2Layers/2 - 1; iLayer++) {  // shift one layer down
    //
    // Place inner absorber layer first
    //
    theta = hcal1DPhi/nHcal1Layers * iLayer;  // add an off set 0.01 rad
    G4double xposShift = rAbsLayerCenter * (cos(theta) - 1.0);
    G4double yposShift = rAbsLayerCenter * sin(theta);
    xPadding = 0.005*RmidAbsLayerX*LayerNum - 2.54*1.9*cm;
    G4ThreeVector absTrans = G4ThreeVector(-RmidAbsLayerX*cos(theta) + xposShift - xPadding, yposShift + yPadding, 0); 
    G4RotationMatrix rotAbsLayer  = G4RotationMatrix();
    rotAbsLayer.rotateY(90*deg);
    rotAbsLayer.rotateX(90*deg);
    rotAbsLayer.rotateY(-90*deg);
    tiltPadding = LayerNum*0.005*rad;
    rotAbsLayer.rotateZ((theta-theta2)*rad + tiltPadding);
    //rotAbsLayer.rotateZ(-(iLayer + nHcal2Layers/2)*0.02*rad + tiltPadding);
    
    G4Transform3D transformAbs = G4Transform3D(rotAbsLayer,absTrans);
    
    new G4PVPlacement(transformAbs,              //rotation,position
		      logicHcal1AbsLayer,        //its logical volume
		      "hcal1AbsLayer",           //its name
		      logicHcalBox,              //its mother  volume
		      false,                     //no boolean operation
		      LayerNum - 1,              //copy number
		      checkOverlaps);            //overlaps checking

    //
    // Place inner hcal scintillator layer
    //

    theta += 0.5*hcal1DPhi/nHcal1Layers;
    //G4cout << "M_PI: " << M_PI << "     theta: " << theta << "    TileAngle: " << theta2 << G4endl;
    //G4double phi = iLayer*dPhi;
    //    G4ThreeVector myTrans = G4ThreeVector(-RmidAbsLayerX*cos(theta) + xposShift,  yposShift, 0); 
    xposShift = rAbsLayerCenter * (cos(theta) - 1.0);
    yposShift = rAbsLayerCenter * sin(theta);
    G4ThreeVector myTrans = G4ThreeVector(-RmidAbsLayerX*cos(theta) + xposShift - xPadding, yposShift + yPadding, 0); 
    //G4cout << " x: " << Rmid*cos(theta) << "    y: " << Rmid*sin(theta) << "   1.0*rad: " << 1.8*rad << G4endl;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateZ((theta-theta2+0.01)*rad + tiltPadding);        // Added a funge increment factor 0.01
    //rotm.rotateZ(-(iLayer + nHcal2Layers/2+1)*0.02*rad + tiltPadding);        // Added a funge increment factor 0.01
    G4Transform3D transform = G4Transform3D(rotm,myTrans);
     
    G4cout << "  iLayer " << iLayer << G4endl;

    new G4PVPlacement(transform,               //rotation,position
                      logicHcal1ScintLayer,    //its logical volume
                      "hcal1ScintLayer",       //its name
                      logicHcalBox,            //its mother volume: logicHcal1Ab
                      false,                   //no boolean operation
                      LayerNum - 1,            //copy number
                      checkOverlaps);          //checking overlaps 
    LayerNum++;
  }

  /*
  theta = 0.;
  theta2 = hcal1TiltAngle;
  G4double Rmid = hcal1RadiusIn + hcal1ScintSizeX / 2.0 + 0.5*cm;  // add extra 0.5cm 
  for (G4int iLayer = 0; iLayer < nHcal1Layers; iLayer++) {  // need to figure out why the index starts at 2 ???
    theta = hcal1DPhi/nHcal1Layers * iLayer - 0.008;  // add an off set 0.008 rad
    //G4cout << "M_PI: " << M_PI << "     theta: " << theta << "    TileAngle: " << theta2 << G4endl;
    //    G4double phi = iLayer*dPhi;
    G4ThreeVector myTrans = G4ThreeVector(Rmid*cos(theta), Rmid*sin(theta), 0); 
    //G4cout << " x: " << Rmid*cos(theta) << "    y: " << Rmid*sin(theta) << "   1.0*rad: " << 1.8*rad << G4endl;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateZ((theta-theta2+0.01)*rad);        // Added a funge increment factor 0.01
    G4Transform3D transform = G4Transform3D(rotm,myTrans);
     
    G4cout << "  iLayer " << iLayer << G4endl;

    new G4PVPlacement(transform,             //rotation,position
                      logicHcal1ScintLayer,           //its logical volume
                      "hcal1ScintLayer",             //its name
                      logicWorld,             //its mother volume: logicHcal1Ab
                      false,                 //no boolean operation
                      iLayer,                 //copy number
                      checkOverlaps);       // checking overlaps 

    theta -= 0.5*hcal1DPhi/nHcal1Layers;  // add an off set 05 rad 
    G4ThreeVector absTrans = G4ThreeVector(Rmid*cos(theta), Rmid*sin(theta), 0); 
    G4RotationMatrix rotAbsLayer  = G4RotationMatrix();
    rotAbsLayer.rotateY(90*deg);
    rotAbsLayer.rotateX(90*deg);
    rotAbsLayer.rotateY(-90*deg);
    rotAbsLayer.rotateZ((theta-theta2)*rad);
    
    G4Transform3D transformAbs = G4Transform3D(rotAbsLayer,absTrans);
    
    new G4PVPlacement(transformAbs,              //rotation,position
		      logicHcal1AbsLayer,        //its logical volume
		      "hcal1AbsLayer",           //its name
		      logicWorld,                //its mother  volume
		      false,                     //no boolean operation
		      iLayer,                    //copy number
		      checkOverlaps);            //overlaps checking
		      }
  
  */							     
  return physiWorld;
}


/*
// here we build up the outer hcal from the sandwiches (steel + scintillator)
int
PHG4HcalPrototypeDetector::ConstructOuterHcal( G4LogicalVolume* logicWorld )
{
  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4VSolid*  outer_sandwich = new G4Trap("Hcal_Outer_Sandwich_Trap", outer_steel_z, outer_steel_x, outer_steel_y_out+ sc_dimension[1]+4*no_overlap, outer_steel_y_in+ sc_dimension[1]+4*no_overlap);
  G4LogicalVolume* outer_sandwich_log =  new G4LogicalVolume(outer_sandwich, Air, "Hcal_Outer_Sandwich", 0, 0, 0);
  G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Blue());
outer_sandwich_log->SetVisAttributes(cemcVisAtt);
  ostringstream outer_hcal_name;
      // hcal_rotm->rotateZ(270 * deg);
      // new G4PVPlacement(hcal_rotm, G4ThreeVector(0, 0, 0), outer_sandwich_log,"hcal_1" , logicWorld, 0, i, overlapcheck);
      // hcal_rotm->rotateZ(90 * deg+i*outer_tilt_angle*rad + 2*deg);
      // new G4PVPlacement(hcal_rotm, G4ThreeVector(0, 0, 0), outer_sandwich_log,"hcal_1" , logicWorld, 0, i, overlapcheck);
      // return 0;
  //  for (int i=0; i<15; i++)
  double yup = -((num_sandwiches-1)*(outer_steel_y_out+sc_dimension[1])+outer_steel_y_out)/2.;
  double xup = 0.;
  for (int i=0; i<15; i++)
    {
      outer_hcal_name.str("");
      outer_hcal_name << "OuterHcal_" << i;
      G4RotationMatrix *hcal_rotm = new G4RotationMatrix;

      hcal_rotm->rotateZ(270 * deg-i*outer_tilt_angle*rad);

      new G4PVPlacement(hcal_rotm, G4ThreeVector(xup, yup, 0), outer_sandwich_log,outer_hcal_name.str().c_str() , logicWorld, 0, i, overlapcheck);
      yup += ((outer_steel_y_out+outer_steel_y_in)/2.+sc_dimension[1]+9*no_overlap)*cos(i*outer_tilt_angle);
      xup -=  ((outer_steel_y_out+outer_steel_y_in)/2.+sc_dimension[1]+9*no_overlap)*sin(i*outer_tilt_angle);
    }
  ConstructSandwichVolume(outer_sandwich_log);

  return 0;
}

// here we put a single tungsten + scintillator sandwich together
int
PHG4HcalPrototypeDetector::ConstructSandwichVolume(G4LogicalVolume* sandwich_log)
{
  vector<G4LogicalVolume*> block_logic;
  G4VSolid*  outer_steel_trd = new G4Trap("Hcal_outer_steel_trd", outer_steel_z, outer_steel_x, outer_steel_y_out, outer_steel_y_in);
  G4LogicalVolume* outer_steel_log =  new G4LogicalVolume(outer_steel_trd, G4Material::GetMaterial("G4_Fe"), "Hcal_Outer_Steel", 0, 0, 0);
  G4VisAttributes* cemcVisAtt = new G4VisAttributes();
  cemcVisAtt->SetVisibility(true);
  cemcVisAtt->SetForceSolid(false);
  cemcVisAtt->SetColour(G4Colour::Magenta());
  outer_steel_log->SetVisAttributes(cemcVisAtt);
  G4VSolid *outer_scintillator = new G4Box("Outer_Scinti_Box", sc_dimension[1] / 2., sc_dimension[0] / 2., sc_dimension[2] / 2.);
  G4LogicalVolume* outer_scinti_log = new G4LogicalVolume(outer_scintillator, G4Material::GetMaterial("G4_POLYSTYRENE"), "Hcal_Outer_Sc", 0, 0, 0);
  G4VisAttributes* cemcVisAtt1 = new G4VisAttributes();
  cemcVisAtt1->SetVisibility(true);
  cemcVisAtt1->SetForceSolid(false);
  cemcVisAtt1->SetColour(G4Colour::Green());
  outer_scinti_log->SetVisAttributes(cemcVisAtt1);
  //  G4RotationMatrix *hcal_rotm = new G4RotationMatrix;
  //  hcal_rotm->rotateZ(270 * deg);
  // get to the bottom of the sandwich volome
  // -((y_out + y_in)/2)/2-sc/2
  // going 1/2 trapezoid up
  // -((y_out + y_in)/2)/2-sc/2 + ((y_out + y_in)/2)/2 = -sc/2

  new G4PVPlacement(NULL, G4ThreeVector( - sc_dimension[1] / 2., 0, 0), outer_steel_log, "HcalSteel", sandwich_log, 0, false, overlapcheck);
  //  new G4PVPlacement(hcal_rotm, G4ThreeVector(0, -((outer_steel_y_out + outer_steel_y_in)/2.)/2. - sc_dimension[1] / 2., 0), outer_steel_log, "HcalSteel", sandwich_log, 0, false, overlapcheck);
  G4RotationMatrix *hcal_rotm_sc = new G4RotationMatrix;
  hcal_rotm_sc->rotateZ((2 * M_PI - outer_tilt_angle)*rad);
  new G4PVPlacement(hcal_rotm_sc, G4ThreeVector(((outer_steel_y_out + outer_steel_y_in)/2.)/2.+no_overlap, 0, 0), outer_scinti_log, "HcalScinti", sandwich_log, 0, false, overlapcheck);
  return 0;
}

*/

void
PHG4HcalPrototypeDetector::CalculateGeometry()
{
  // tilt angle outer steel
  // outer_tilt_angle = atan((outer_steel_y_out-outer_steel_y_in)/outer_plate_x);
  // cout << "outer tilt angle: " << outer_tilt_angle*180./M_PI << endl;
    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PHG4HcalPrototypeDetector::SetMaterial(G4String materialChoice)
{
  /*
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (pttoMaterial) {
    world_mat = pttoMaterial;
    logicWorld->SetMaterial(world_mat);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void PHG4HcalPrototypeDetector::UpdateGeometry()
{
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
}

// The following code is copied from the sPHENIX G4 simulation
void PHG4HcalPrototypeDetector::SetTiltViaNcross(const int ncross)
{
  G4double sign = 1;
  if (ncross < 0)
    {
      sign = -1;
    }
  G4int ncr = fabs(ncross);

  // Determine title angle for the outer hcal
  //
  G4double cSide = hcal2RadiusIn + hcal2ScintSizeX/2;
  G4double bSide = hcal2RadiusIn + hcal2ScintSizeX;

  G4double alpha = 0;
  if (ncr > 1)
    {
      alpha = (360. / nHcal2Layers * M_PI / 180.) * (ncr - 1) / 2.0;
    }
  else
    {
      alpha = (360. / nHcal2Layers * M_PI / 180.) / 2.;
    }

  G4double sinbSide = sin(alpha) * bSide / (sqrt(bSide * bSide + cSide * cSide - 2 * bSide * cSide * cos(alpha)));
  G4double beta = asin(sinbSide);  // This is the slat angle

  hcal2TiltAngle = beta * sign;

  // Determine title angle for the inner hcal
  //
  cSide = hcal1RadiusIn + hcal1ScintSizeX/2;
  bSide = hcal1RadiusIn + hcal1ScintSizeX;

  /* This part is the same both for the inner and the outer hcal section
  G4doube alpha = 0;
  if (ncr > 1)
    {
      alpha = (360. / nScint360 * M_PI / 180.) * (ncr - 1);
    }
  else
    {
      alpha = (360. / nScint360 * M_PI / 180.) / 2.;
    }
  */

  sinbSide = sin(alpha) * bSide / (sqrt(bSide * bSide + cSide * cSide - 2.0 * bSide * cSide * cos(alpha)));
  beta = asin(sinbSide);  // This is the slat angle

  hcal1TiltAngle = beta * sign;

  G4cout << " alpha : " << alpha << G4endl;
  G4cout <<  " SetTitlViaNCross(" << ncross << ") setting the outer hcal slat tilt angle to : " << hcal2TiltAngle << " radian" << G4endl;
  G4cout <<  " SetTitlViaNCross(" << ncross << ") setting the inner hcal slat tilt angle to : " << hcal1TiltAngle << " radian" << G4endl;
  return;

}

// Detector construction messengers for setting plate angles
// 
void PHG4HcalPrototypeDetector::SetOuterHcalDPhi(G4double dphi)
{
  hcal2DPhi = dphi;
  G4cout << "In SetOuterHcalDPhi: " << hcal2DPhi << " is set!!! " << G4endl;
}

void PHG4HcalPrototypeDetector::SetOuterPlateTiltAngle(G4double dtheta)
{
  hcal2TiltAngle = dtheta;
  G4cout << "In SetOuterPlateTiltAngle: " << hcal2TiltAngle << " is set!!! " << G4endl;
}

void PHG4HcalPrototypeDetector::SetInnerHcalDPhi(G4double dphi)
{
  hcal1DPhi = dphi;
  G4cout << "In SetInnerHcalDPhi: " << hcal1DPhi << " is set!!! " << G4endl;
}

void PHG4HcalPrototypeDetector::SetInnerPlateTiltAngle(G4double dtheta)
{
  hcal1TiltAngle = dtheta;
  G4cout << "In SetInnerPlateTiltAngle: " << hcal1TiltAngle << " is set!!! " << G4endl;
}


