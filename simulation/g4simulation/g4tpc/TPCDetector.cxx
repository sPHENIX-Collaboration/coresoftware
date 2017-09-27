#include <iostream>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4UserLimits.hh>

#include <phool/PHCompositeNode.h>
#include <TPCConstants.h>

#include <TString.h>

#include "TPCDetector.h"

using namespace TPCConstants;

TPCDetector::TPCDetector(PHCompositeNode *node) : 
  PHG4Detector(node),
  fVerbosity(0),
  fTPC(NULL),
  fLimits( new G4UserLimits(0.3*cm) )
{
}
//=====
void TPCDetector::BuildFiducial(G4LogicalVolume* tpcWorld)
{
 
  G4Material *mat_gas = G4Material::GetMaterial("sPHENIX_tpc_gas");

  G4double minR = kGasInnerRadius*cm;
  G4double maxR = kGasOuterRadius*cm;
  G4double hz = kGasHalfOfLength*cm;

  G4Tubs *shp = new G4Tubs( "TPC_GAS_SL", minR, maxR, hz, 0*deg, 360*deg );
  G4LogicalVolume *lv = new G4LogicalVolume( shp, mat_gas, "TPC_GAS_LL" );
  lv->SetUserLimits(fLimits);
  fTPC = new G4PVPlacement(0, G4ThreeVector(0,0,0), lv, "TPC_GAS_PV", tpcWorld, 0, false, overlapcheck );
  /*
  const int lyrs = 1;
  G4double deltar = (maxR-minR)/lyrs;
  G4Tubs *shp[lyrs];
  G4LogicalVolume *lv[lyrs];
  for(int i=0; i!=lyrs; ++i) {
    G4double rmin = minR + i*deltar;
    G4double rmax = minR + (i+1)*deltar;
    shp[i] = new G4Tubs( Form("TPC_GAS_SL%d",i), rmin, rmax, hz, 0*deg, 360*deg );
    lv[i] = new G4LogicalVolume( shp[i], mat_gas, Form("TPC_GAS_LL%d",i) );
    lv[i]->SetUserLimits(fLimits);
    fTPC = new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[i], Form("TPC_GAS_PV%d",i), tpcWorld, 0, false, overlapcheck );
  }
  */
}
//=====
void TPCDetector::BuildCage(G4LogicalVolume* tpcWorld)
{
   
  G4Material *mat_cu = G4Material::GetMaterial("G4_Cu");
  G4Material *mat_kap = G4Material::GetMaterial("G4_KAPTON");
  G4Material *mat_fr4 = G4Material::GetMaterial("G4_G10");
  G4Material *mat_hco = G4Material::GetMaterial("G4_MYLAR");

  G4double rads[7];
  G4double hz = kGasHalfOfLength*cm;
  G4Tubs *shp[6];
  G4LogicalVolume *lv[6];

  //-- INNER
  rads[0] = (kGasInnerRadius-1.17)*cm;
  rads[1] = rads[0] +  0.05*mm; // 50 um Cu
  rads[2] = rads[1] +  0.05*mm; // 50 um FR4
  rads[3] = rads[2] + 10.00*mm; // 1 cm HONEYCOMB
  rads[4] = rads[3] +  0.05*mm; // 50 um FR4
  rads[5] = rads[4] +  1.50*mm; // 1.5 mm KAPTON
  rads[6] = rads[5] +  0.05*mm; // 50 um FR4

  shp[0] = new G4Tubs( "TPC_IC_SL0", rads[0], rads[1], hz, 0*deg, 360*deg );
  shp[1] = new G4Tubs( "TPC_IC_SL1", rads[1], rads[2], hz, 0*deg, 360*deg );
  shp[2] = new G4Tubs( "TPC_IC_SL2", rads[2], rads[3], hz, 0*deg, 360*deg );
  shp[3] = new G4Tubs( "TPC_IC_SL3", rads[3], rads[4], hz, 0*deg, 360*deg );
  shp[4] = new G4Tubs( "TPC_IC_SL4", rads[4], rads[5], hz, 0*deg, 360*deg );
  shp[5] = new G4Tubs( "TPC_IC_SL5", rads[5], rads[6], hz, 0*deg, 360*deg );

  lv[0] = new G4LogicalVolume( shp[0], mat_cu,  "TPC_IC_LL0");
  lv[1] = new G4LogicalVolume( shp[1], mat_fr4, "TPC_IC_LL1");
  lv[2] = new G4LogicalVolume( shp[2], mat_hco, "TPC_IC_LL2");
  lv[3] = new G4LogicalVolume( shp[3], mat_fr4, "TPC_IC_LL3");
  lv[4] = new G4LogicalVolume( shp[4], mat_kap, "TPC_IC_LL4");
  lv[5] = new G4LogicalVolume( shp[5], mat_fr4, "TPC_IC_LL5");

  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[0], "TPC_IC_PL0", tpcWorld, 0, false, overlapcheck );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[1], "TPC_IC_PL1", tpcWorld, 0, false, overlapcheck );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[2], "TPC_IC_PL2", tpcWorld, 0, false, overlapcheck );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[3], "TPC_IC_PL3", tpcWorld, 0, false, overlapcheck );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[4], "TPC_IC_PL4", tpcWorld, 0, false, overlapcheck );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[5], "TPC_IC_PL5", tpcWorld, 0, false, overlapcheck );

  //-- OUTER

  rads[0] = kGasOuterRadius*cm;
  rads[1] = rads[0] + 0.050*mm; // 50 um Cu
  rads[2] = rads[1] + 0.050*mm; // 50 um FR4
  rads[3] = rads[2] + 10.00*mm; // 1 cm HONEYCOMB
  rads[4] = rads[3] + 0.050*mm; // 50 um FR4
  rads[5] = rads[4] + 1.500*mm; // 1.5 mm KAPTON
  rads[6] = rads[5] + 0.050*mm; // 50 um PCB

  shp[0] = new G4Tubs( "TPC_OC_SL0", rads[0], rads[1], hz, 0*deg, 360*deg );
  shp[1] = new G4Tubs( "TPC_OC_SL1", rads[1], rads[2], hz, 0*deg, 360*deg );
  shp[2] = new G4Tubs( "TPC_OC_SL2", rads[2], rads[3], hz, 0*deg, 360*deg );
  shp[3] = new G4Tubs( "TPC_OC_SL3", rads[3], rads[4], hz, 0*deg, 360*deg );
  shp[4] = new G4Tubs( "TPC_OC_SL4", rads[4], rads[5], hz, 0*deg, 360*deg );
  shp[5] = new G4Tubs( "TPC_OC_SL5", rads[5], rads[6], hz, 0*deg, 360*deg );

  lv[0] = new G4LogicalVolume( shp[0], mat_cu,  "TPC_OC_LL0");
  lv[1] = new G4LogicalVolume( shp[1], mat_fr4, "TPC_OC_LL1");
  lv[2] = new G4LogicalVolume( shp[2], mat_hco, "TPC_OC_LL2");
  lv[3] = new G4LogicalVolume( shp[3], mat_fr4, "TPC_OC_LL3");
  lv[4] = new G4LogicalVolume( shp[4], mat_kap, "TPC_OC_LL4");
  lv[5] = new G4LogicalVolume( shp[5], mat_fr4, "TPC_OC_LL5");

  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[0], "TPC_OC_PL0", tpcWorld, 0, false, overlapcheck );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[1], "TPC_OC_PL1", tpcWorld, 0, false, overlapcheck );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[2], "TPC_OC_PL2", tpcWorld, 0, false, overlapcheck );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[3], "TPC_OC_PL3", tpcWorld, 0, false, overlapcheck );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[4], "TPC_OC_PL4", tpcWorld, 0, false, overlapcheck );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), lv[5], "TPC_OC_PL5", tpcWorld, 0, false, overlapcheck );
}
//=====
void TPCDetector::Construct(G4LogicalVolume* world)
{
  if(fVerbosity>0) std::cout << "TPCDetector::Construct" << std::endl;

  G4Material *mat_air = G4Material::GetMaterial("G4_AIR") ;
  G4double minR = 19.5*cm;
  G4double maxR = 78.5*cm;
  G4double hz = 110*cm;
  G4Tubs *tpc = new G4Tubs( "TPC_SV", minR, maxR, hz, 0*deg, 360*deg );
  G4LogicalVolume *tpc_world = new G4LogicalVolume( tpc, mat_air, "TPC_LV" );
  //BuildCage( tpc_world );
  BuildFiducial( tpc_world );
  new G4PVPlacement(0, G4ThreeVector(0,0,0), tpc_world, "TPC_PV", world, 0, false, overlapcheck );
}
//=====
bool TPCDetector::IsThisActive(G4VPhysicalVolume *test)
{
  if(test==fTPC) return true;
  return false;
}
//=====
void TPCDetector::SetMaxStep(G4double maxStep)
{
  if ((fLimits)&&(maxStep>0.)) fLimits->SetMaxAllowedStep(maxStep);
}
