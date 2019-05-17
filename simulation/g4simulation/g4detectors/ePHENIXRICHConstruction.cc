// $$Id: ePHENIXRICHConstruction.cc,v 1.7 2014/05/01 19:02:45 phnxbld Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov> and Nils Feege <nils.feege@stonybrook.edu>
 * \version $$Revision: 1.7 $$
 * \date $$Date: 2014/05/01 19:02:45 $$
 */

#include "ePHENIXRICHConstruction.h"
#include "PHG4RICHDisplayAction.h"
#include "PHG4RICHSubsystem.h"

#include <Geant4/G4Cons.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4LogicalSkinSurface.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4OpticalSurface.hh>
#include <Geant4/G4Orb.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Sphere.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Tubs.hh>

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;
using namespace ePHENIXRICH;

ePHENIXRICHConstruction::ePHENIXRICHConstruction(PHG4RICHSubsystem *subsys)
  : m_DisplayAction(dynamic_cast<PHG4RICHDisplayAction*>(subsys->GetDisplayAction())),
overlapcheck_rich(false)

{
}
ePHENIXRICHConstruction::ePHENIXRICHConstruction(PHG4RICHSubsystem *subsys, const RICH_Geometry &g)
  :  
  geom(g),
m_DisplayAction(dynamic_cast<PHG4RICHDisplayAction*>(subsys->GetDisplayAction())),
  overlapcheck_rich(false)
{
}

G4LogicalVolume *ePHENIXRICHConstruction::RegisterLogicalVolume(G4LogicalVolume *v)
{
  if (!v)
  {
    G4cout << "ePHENIXRICHConstruction::RegisterVolume - Error - invalid volume!"
           << G4endl;
    return v;
  }
  if (map_log_vol.find(v->GetName()) != map_log_vol.end())
  {
    G4cout << "ePHENIXRICHConstruction::RegisterVolume - Warning - replacing "
           << v->GetName() << G4endl;
  }

  map_log_vol[v->GetName()] = v;

  return v;
}

G4PVPlacement *ePHENIXRICHConstruction::RegisterPhysicalVolume(G4PVPlacement *v)
{
  if (!v)
  {
    G4cout << "ePHENIXRICHConstruction::RegisterPhysicalVolume - Error - invalid volume!"
           << G4endl;
    return v;
  }

  phy_vol_idx_t id(v->GetName(), v->GetCopyNo());

  if (map_phy_vol.find(id) != map_phy_vol.end())
  {
    G4cout << "ePHENIXRICHConstruction::RegisterPhysicalVolume - Warning - replacing "
           << v->GetName() << "[" << v->GetCopyNo() << "]" << G4endl;
  }

  map_phy_vol[id] = v;
  return v;
}

G4LogicalVolume *
ePHENIXRICHConstruction::Construct_RICH(G4LogicalVolume *WorldLog)
{
  // -- Logical volume:
  G4VSolid *RICHOutSphereBoundary = new G4Orb("RICHOutSphereBoundary",
                                              geom.get_R_max());
  G4VSolid *RICHOutSphereBoundary_place = new G4DisplacedSolid("RICHOutSphereBoundary_place", RICHOutSphereBoundary, 0,
                                                               G4ThreeVector(0, geom.get_R_shift(), geom.get_z_shift()));

  G4VSolid *RICHInnerSphereBoundary = new G4Orb("RICHInnerSphereBoundary",
                                                geom.get_R_frontwindow());
  G4VSolid *RICHInnerSphereBoundary_place = new G4DisplacedSolid("RICHInnerSphereBoundary_place",
                                                                 RICHInnerSphereBoundary,
                                                                 0,
                                                                 G4ThreeVector(0,
                                                                               geom.get_R_shift() * geom.get_frontwindow_DisplaceRatio(),
                                                                               geom.get_z_shift() * geom.get_frontwindow_DisplaceRatio()));

  G4VSolid *RICHConeBoundary = new G4Cons("RICHConeBoundary",  //
                                          geom.get_R_beam_pipe(),
                                          geom.get_z_shift() / 2 * std::tan(2 * std::atan(std::exp(-geom.get_min_eta()))),  //            G4double pRmin1, G4double pRmax1,
                                          geom.get_R_beam_pipe(),
                                          geom.get_cone_size_z() * std::tan(2 * std::atan(std::exp(-geom.get_min_eta()))),  //            G4double pRmin2, G4double pRmax2,
                                          (geom.get_cone_size_z() - (geom.get_z_shift() / 2)) / 2,                          //            G4double pDz,
                                          0, 2 * pi                                                                         //            G4double pSPhi, G4double pDPhi
                                          );
  G4VSolid *RICHConeBoundary_place = new G4DisplacedSolid("RICHConeBoundary_place",
                                                          RICHConeBoundary,
                                                          0,
                                                          G4ThreeVector(
                                                              0,
                                                              0,
                                                              (geom.get_cone_size_z() - (geom.get_z_shift() / 2)) / 2 + (geom.get_z_shift() / 2)));

  G4VSolid *RICHSecBoundary = new G4Tubs("RICHSecBoundary",                       //
                                         geom.get_R_beam_pipe(),                  //            G4double pRMin,
                                         geom.get_cone_size_z(),                  //            G4double pRMax,
                                         geom.get_cone_size_z(),                  //            G4double pDz,
                                         pi / 2 - pi / geom.get_N_RICH_Sector(),  //            G4double pSPhi,
                                         2 * pi / geom.get_N_RICH_Sector()        //            G4double pDPhi
                                         );

  G4VSolid *RICHSphereBoundary = new G4SubtractionSolid("RICHSphereBoundary",
                                                        RICHOutSphereBoundary_place, RICHInnerSphereBoundary_place);
  //      G4VSolid *RICHSecBox_ConeSphere = new G4IntersectionSolid(
  //              "RICHSecBox_ConeSphere", RICHSphereBoundary,
  //              RICHConeBoundary_place);
  G4VSolid *RICHSecBox_ConeSec = new G4IntersectionSolid("RICHSecBox_ConeSec",
                                                         RICHSecBoundary, RICHConeBoundary_place);

  // RICH sector
  G4VSolid *RICHSecBox = new G4IntersectionSolid("RICHSecBox",
                                                 RICHSecBox_ConeSec, RICHSphereBoundary);

  G4LogicalVolume *RICHSecLog = new G4LogicalVolume(RICHSecBox,
                                                    G4Material::GetMaterial(geom.get_RICH_gas_mat()), "RICHSecLogical",
                                                    0, 0, 0);
  RegisterLogicalVolume(RICHSecLog);

  for (G4int sec = 0; sec < geom.get_N_RICH_Sector(); sec++)
  {
    G4RotateZ3D sec_rot(2 * pi / geom.get_N_RICH_Sector() * sec);

    // Insert sector volumes into sector set
    sector_vec.insert(RegisterPhysicalVolume(new G4PVPlacement(sec_rot, RICHSecLog, "RICHSecPhysical", WorldLog,
                                                               false, sec, overlapcheck_rich)));
  }
  //      CLHEP::HepRotationZ sec_rot(0);
  //      G4RotationMatrix g4_sec_rot(sec_rot);
  //      new G4PVPlacement(transform1, "RICHSecPhysical", RICHSecLog, WorldPhys,
  //              false, 0);

  //RICH mirror
  // LHCb website: 3mm thick Be base + 0.3mm glass surface layer coated with Al.
  // LHCb TDR: The mirrors are made of polished 6mm-thick glass coated by vacuum
  //           deposition with 900 nm of aluminium and overcoated with 200 nm of quartz.
  G4VSolid *RICHMirrorSphereBoundary = new G4Sphere("RICHMirrorSphereBoundary",  //
                                                    geom.get_R_mirror_ref(),
                                                    geom.get_R_mirror_ref() + geom.get_dR_mirror(),  //            G4double pRmin, G4double pRmax,
                                                    0, 2 * pi,                                       //            G4double pSPhi, G4double pDPhi,
                                                    0, pi                                            //            G4double pSTheta, G4double pDTheta
                                                    );
  G4VSolid *RICHMirrorSphereBoundary_place = new G4DisplacedSolid("RICHMirrorSphereBoundary_place", RICHMirrorSphereBoundary, 0,
                                                                  G4ThreeVector(0, geom.get_R_shift(), geom.get_z_shift()));
  G4VSolid *RICHMirror = new G4IntersectionSolid("RICHMirror",
                                                 RICHSecBox_ConeSec, RICHMirrorSphereBoundary_place);
  G4LogicalVolume *RICHMirrorLog = new G4LogicalVolume(RICHMirror,
                                                       G4Material::GetMaterial(geom.get_RICH_Mirror_mat()),
                                                       "RICHMirrorLog");
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(), RICHMirrorLog, "RICHMirrorPhysical",
                                           RICHSecLog, false, 0, overlapcheck_rich));
  RegisterLogicalVolume(RICHMirrorLog);

  // RICH mirror Optical Surface
  //
  G4LogicalSkinSurface *RICHMirrorSurface = new G4LogicalSkinSurface("RICHMirrorSurface",
                                                                     RICHMirrorLog,
                                                                     geom.get_RICH_Mirror_OpticalSurface());
  RICHMirrorSurface->GetName();  // only to avoid error message for unused object

  //RICH gas vessel window - back
  //LHCb TDR: The frame will be sealed to contain the C4F10 gas radiator. The vacuum chamber acts as part of the boundary to the gas volume.
  //                    Kapton foils of 150um thickness and 400 mm diameter will be glued to flanges on the vacuum chamber
  G4VSolid *RICHBackWindowSphereBoundary = new G4Sphere("RICHBackWindowSphereBoundary",  //
                                                        geom.get_R_mirror_ref() + geom.get_dR_mirror() + geom.get_dR_mirror_spt(),
                                                        geom.get_R_mirror_ref() + geom.get_dR_mirror() + geom.get_dR_mirror_spt() + geom.get_dR_backwindow(),  //            G4double pRmin, G4double pRmax,
                                                        0, 2 * pi,                                                                                             //            G4double pSPhi, G4double pDPhi,
                                                        0, pi                                                                                                  //            G4double pSTheta, G4double pDTheta
                                                        );
  G4VSolid *RICHBackWindowSphereBoundary_place = new G4DisplacedSolid("RICHBackWindowSphereBoundary_place", RICHBackWindowSphereBoundary,
                                                                      0, G4ThreeVector(0, geom.get_R_shift(), geom.get_z_shift()));
  G4VSolid *RICHBackWindow = new G4IntersectionSolid("RICHBackWindow",
                                                     RICHSecBox_ConeSec, RICHBackWindowSphereBoundary_place);
  G4LogicalVolume *RICHBackWindowLog = new G4LogicalVolume(RICHBackWindow,
                                                           G4Material::GetMaterial(geom.get_RICH_Gas_Window_mat()),
                                                           "RICHBackWindowLog");
  RegisterLogicalVolume(RICHBackWindowLog);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(), RICHBackWindowLog,
                                           "RICHBackWindowPhysical", RICHSecLog, false, 0, overlapcheck_rich));

  //RICH gas vessel window - front
  //LHCb TDR: The frame will be sealed to contain the C4F10 gas radiator. The vacuum chamber acts as part of the boundary to the gas volume.
  //                    Kapton foils of 150um thickness and 400 mm diameter will be glued to flanges on the vacuum chamber
  G4VSolid *RICHFrontWindowSphereBoundary = new G4Sphere("RICHFrontWindowSphereBoundary",  //
                                                         geom.get_R_frontwindow(),
                                                         geom.get_R_frontwindow() + geom.get_dR_frontwindow(),  //            G4double pRmin, G4double pRmax,
                                                         0, 2 * pi,                                             //            G4double pSPhi, G4double pDPhi,
                                                         0, pi                                                  //            G4double pSTheta, G4double pDTheta
                                                         );
  G4VSolid *RICHFrontWindowSphereBoundary_place = new G4DisplacedSolid("RICHFrontWindowSphereBoundary_place",
                                                                       RICHFrontWindowSphereBoundary,
                                                                       0,
                                                                       G4ThreeVector(0,
                                                                                     geom.get_R_shift() * geom.get_frontwindow_DisplaceRatio(),
                                                                                     geom.get_z_shift() * geom.get_frontwindow_DisplaceRatio()));
  G4VSolid *RICHFrontWindow = new G4IntersectionSolid("RICHFrontWindow",
                                                      RICHSecBox_ConeSec, RICHFrontWindowSphereBoundary_place);
  G4LogicalVolume *RICHFrontWindowLog = new G4LogicalVolume(RICHFrontWindow,
                                                            G4Material::GetMaterial(geom.get_RICH_Gas_Window_mat()),
                                                            "RICHFrontWindowLog");
  RegisterLogicalVolume(RICHFrontWindowLog);
  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(), RICHFrontWindowLog,
                                           "RICHFrontWindowPhysical", RICHSecLog, false, 0, overlapcheck_rich));

  // photon detector - HBD
  G4LogicalVolume *RICHHBDLog = Construct_HBD(RICHSecLog);

  GetDisplayAction()->AddVolume(RICHSecLog,"Sector");

  GetDisplayAction()->AddVolume(RICHMirrorLog,"Mirror");

  GetDisplayAction()->AddVolume(RICHBackWindowLog,"Window");
  GetDisplayAction()->AddVolume(RICHFrontWindowLog,"Window");

  GetDisplayAction()->AddVolume(RICHHBDLog,"HBD");

  G4cout << "ePHENIXRICHConstruction::Construct_RICH - " << map_log_vol.size()
         << " logical volume constructed" << G4endl;
  G4cout << "ePHENIXRICHConstruction::Construct_RICH - " << map_phy_vol.size()
         << " physical volume constructed" << G4endl;

  return WorldLog;
}

G4LogicalVolume *
ePHENIXRICHConstruction::Construct_HBD(G4LogicalVolume *RICHSecLog)
{
  const double HBD_thickness = geom.get_HBD_thickness();
  const int n_GEM_layers = geom.get_n_GEM_layers();

  assert(
      HBD_thickness < geom.get_dR_frontwindow_shrink() - geom.get_dR_frontwindow());
  // depth check

  G4VSolid *RICHHBDBox = new G4Cons("RICHHBDBox",                                                       //
                                    0, geom.get_RZ_Seg1_HBD() + geom.get_RZ_Seg2_HBD(),                 //            G4double pRmin1, G4double pRmax1,
                                    0, geom.get_RZ_Seg1_HBD() + geom.get_RZ_Seg2_HBD(),                 //            G4double pRmin2, G4double pRmax2,
                                    HBD_thickness / 2,                                                  //            G4double pDz,
                                    -geom.get_half_angle_HBD() + pi / 2, 2 * geom.get_half_angle_HBD()  //            G4double pSPhi, G4double pDPhi
                                    );

  G4LogicalVolume *RICHHBDLog = new G4LogicalVolume(RICHHBDBox,
                                                    G4Material::GetMaterial(geom.get_RICH_gas_mat()), "RICHHBDLog");
  RegisterLogicalVolume(RICHHBDLog);

  G4Transform3D transform1 = G4Translate3D(0, geom.get_R_Tip_HBD(),
                                           geom.get_Z_Tip_HBD()) *
                             G4RotateX3D(-geom.get_Rotation_HBD()) * G4RotateY3D(pi) * G4TranslateZ3D(HBD_thickness / 2);
  //  G4Transform3D transform1 = G4TranslateZ3D(HBD_thickness / 2);

  RegisterPhysicalVolume(new G4PVPlacement(transform1, RICHHBDLog, "RICHHBDPhysical", RICHSecLog,
                                           false, 0, overlapcheck_rich));

  // Internal HBD structure
  // From doi:10.1016/j.nima.2011.04.015
  // Component Material X0 (cm) Thickness (cm) Area (%) Rad. Length (%)
  double current_z = -HBD_thickness / 2;

  //  Mesh SS 1.67 0.003 11.5 0.021
  // exclude for now to let all potons pass to the top GEM foil
  //G4LogicalVolume* HBDMeshLog = Construct_HBD_Layers(RICHHBDLog, "Mesh", "Steel", current_z,
  //                                                   0.003 * cm * 11.5e-2);
  current_z += 0.003 * cm;

  // photon detector - photocathode is CsI coating on top layer (surface)
  //
  G4LogicalSkinSurface *RICHPhotocathodeSurface;

  //  GEM frames FR4 17.1 0.15x4 6.5 0.228 ; 6.5 is the percentage of the area covered by material.
  // exclude for now to let all potons pass to the top GEM foil
  //Construct_HBD_Layers(RICHHBDLog, "Frame0", "G10", current_z,
  //                     0.15 * cm * 6.5e-2);
  current_z += 0.15 * cm;

  for (int gem = 1; gem <= n_GEM_layers; gem++)
  {
    stringstream sid;
    sid << gem;

    //  GEM Copper 1.43 0.0005x6 64 0.134
    G4LogicalVolume *GemTopLog = Construct_HBD_Layers(RICHHBDLog,
                                                      G4String("GEMFrontCu") + G4String(sid.str()), "G4_Cu",
                                                      current_z, 0.0005 * cm * 64e-2);
    current_z += 0.0005 * cm;

    // top GEM has CsI coating for photon detection
    if (gem == 1)
    {
      RICHPhotocathodeSurface = new G4LogicalSkinSurface("RICHPhotocathodeSurface", GemTopLog, geom.get_RICH_Photocathode_OpticalSurface());
      RICHPhotocathodeSurface->GetName();  // avoid error due to unused object
    }

    //  GEM Kapton 28.6 0.005x3 64 0.034
    Construct_HBD_Layers(RICHHBDLog,
                         G4String("GEMKapton") + G4String(sid.str()), "G4_KAPTON",
                         current_z, 0.005 * cm * 64e-2);
    current_z += 0.005 * cm;

    //  GEM Copper 1.43 0.0005x6 64 0.134
    Construct_HBD_Layers(RICHHBDLog,
                         G4String("GEMBackCu") + G4String(sid.str()), "G4_Cu", current_z,
                         0.0005 * cm * 64e-2);
    current_z += 0.0005 * cm;

    //  GEM frames FR4 17.1 0.15x4 6.5 0.228
    Construct_HBD_Layers(RICHHBDLog,
                         G4String("Frame") + G4String(sid.str()), "G10", current_z,
                         0.15 * cm * 6.5e-2);
    current_z += 0.15 * cm;
  }

  //  PCB Kapton 28.6 0.005 100 0.017
  Construct_HBD_Layers(RICHHBDLog, G4String("PCBKapton"), "G4_KAPTON",
                       current_z, 0.005 * cm * 100e-2);
  current_z += 0.005 * cm;

  //  PCB Copper 1.43 0.0005 80 0.028
  Construct_HBD_Layers(RICHHBDLog, G4String("PCBCu"), "G4_Cu", current_z,
                       0.0005 * cm * 80e-2);
  current_z += 0.0005 * cm;

  //  Facesheet FR4 17.1 0.025x2 100 0.292
  Construct_HBD_Layers(RICHHBDLog, "Facesheet", "G10", current_z,
                       0.025 * 2 * cm * 100e-2);
  current_z += 0.025 * 2 * cm;

  //  Panel core Honeycomb 8170 1.905 100 0.023 <- very thin-X0 stuff, ignore

  //  Total vessel 0.82
  //  Readout
  //  Readout board FR4/copper 17.1/1.43 0.05/0.001 100 0.367
  Construct_HBD_Layers(RICHHBDLog, G4String("ReadoutFR4"), "G10", current_z,
                       0.05 * cm * 100e-2);
  current_z += 0.05 * cm;
  Construct_HBD_Layers(RICHHBDLog, G4String("ReadoutCu"), "G4_Cu", current_z,
                       0.001 * cm * 100e-2);
  current_z += 0.001 * cm;

  //  Preamps + sockets Copper 1.43 0.0005 100 0.66
  //  Total readout 1.03
  Construct_HBD_Layers(RICHHBDLog, G4String("SocketsCu"), "G4_Cu", current_z,
                       0.0005 * cm * 100e-2);
  current_z += 0.0005 * cm;

  assert(current_z < HBD_thickness / 2);
  //boundary check

  return RICHHBDLog;
}

G4LogicalVolume *
ePHENIXRICHConstruction::Construct_HBD_Layers(G4LogicalVolume *RICHHBDLog,
                                              const G4String name, const G4String material, const double start_z,
                                              const double thickness)
{
  G4VSolid *box = new G4Cons(G4String("RICHHBD") + name + G4String("Box"),                       //
                             0, geom.get_RZ_Seg1_HBD() + geom.get_RZ_Seg2_HBD(),                 //            G4double pRmin1, G4double pRmax1,
                             0, geom.get_RZ_Seg1_HBD() + geom.get_RZ_Seg2_HBD(),                 //            G4double pRmin2, G4double pRmax2,
                             thickness / 2,                                                      //            G4double pDz,
                             -geom.get_half_angle_HBD() + pi / 2, 2 * geom.get_half_angle_HBD()  //            G4double pSPhi, G4double pDPhi
                             );

  G4LogicalVolume *Log = new G4LogicalVolume(box,
                                             G4Material::GetMaterial(material),  //
                                             G4String("RICHHBD") + name + G4String("Log"));
  RegisterLogicalVolume(Log);

  RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, 0, start_z + thickness / 2), Log,
                                           G4String("RICHHBD") + name + G4String("Physical"), RICHHBDLog,
                                           false, 0, overlapcheck_rich));

  return Log;
}

void RICH_Geometry::SetDefault()
{
  N_RICH_Sector = 8;
  min_eta = 1;
  R_beam_pipe = 3 * cm;
  z_shift = 100 * cm;
  R_shift = 40 * cm;
  frontwindow_DisplaceRatio = .85;  // Displace R,Z and radius simultainously
  dR_frontwindow_shrink = 2 * cm;
  R_mirror_ref = 200 * cm;
  dR_mirror = 0.6 * cm;
  dR_mirror_spt = 1.4 * cm;
  dR_backwindow = 150 * um;
  dR_frontwindow = 150 * um;

  HBD_thickness = 1.5 * cm;
  n_GEM_layers = 5;

  RICH_gas_mat = "CF4";
  RICH_Mirror_mat = "G4_Pyrex_Glass";
  RICH_Gas_Window_mat = "G4_KAPTON";
}

void RICH_Geometry::CreateOpticalSurfaces()
{
  // Mirror
  //
  delete RICH_Mirror_OpticalSurface;
  RICH_Mirror_OpticalSurface = new G4OpticalSurface("RICHMirrorSurfaceOptical");

  const G4int NUM = 2;

  G4double pp[NUM] = {2.038 * eV, 4.144 * eV};
  G4double rindex[NUM] = {1.4, 1.4};
  G4double reflectivity[NUM] = {1.0, 1.0};
  G4double efficiency[NUM] = {0.0, 0.0};

  G4MaterialPropertiesTable *RICH_Mirror_OpticalSurface_SMPT = new G4MaterialPropertiesTable();

  RICH_Mirror_OpticalSurface_SMPT->AddProperty("RINDEX", pp, rindex, NUM);
  RICH_Mirror_OpticalSurface_SMPT->AddProperty("REFLECTIVITY", pp, reflectivity, NUM);
  RICH_Mirror_OpticalSurface_SMPT->AddProperty("EFFICIENCY", pp, efficiency, NUM);

  RICH_Mirror_OpticalSurface->SetType(dielectric_metal);
  RICH_Mirror_OpticalSurface->SetModel(glisur);
  RICH_Mirror_OpticalSurface->SetFinish(polished);

  // Photocathode
  //
  RICH_Photocathode_OpticalSurface = new G4OpticalSurface("RICH_Photocathode_OpticalSurface", glisur, polished, dielectric_metal);

  G4double photocath_EPHOTON[2] = {1., 1000.};
  G4double photocath_EFF[2] = {1., 1.};
  G4double photocath_REFL[2] = {0., 0.};

  G4MaterialPropertiesTable *RICH_Photocathode_OpticalSurface_MPT = new G4MaterialPropertiesTable();
  RICH_Photocathode_OpticalSurface_MPT->AddProperty("EFFICIENCY", photocath_EPHOTON, photocath_EFF, 2);
  RICH_Photocathode_OpticalSurface_MPT->AddProperty("REFLECTIVITY", photocath_EPHOTON, photocath_REFL, 2);
  RICH_Photocathode_OpticalSurface->SetMaterialPropertiesTable(RICH_Photocathode_OpticalSurface_MPT);

  // done
  //
  return;
}

double RICH_Geometry::get_R_frontwindow() const
{
  return (R_mirror_ref / 2 + sqrt(z_shift * z_shift + R_shift * R_shift) * (1 - frontwindow_DisplaceRatio)) - dR_frontwindow_shrink;
}

double RICH_Geometry::get_half_angle_HBD() const
{
  return atan(
      tan(pi / N_RICH_Sector) / sqrt(1 + R_shift * R_shift / z_shift / z_shift));
}

double RICH_Geometry::get_RZ_Seg1_HBD() const
{
  return R_shift * (R_mirror_ref + 2 * sqrt(z_shift * z_shift + R_shift * R_shift)) / 4 / z_shift;
}

double RICH_Geometry::get_RZ_Seg2_HBD() const
{
  return ((R_mirror_ref + 2 * sqrt(pow(z_shift, 2) + pow(R_shift, 2))) * tan(2 * atan(exp(-min_eta)) - atan(R_shift / z_shift))) / 4.;
}

double RICH_Geometry::get_R_Tip_HBD() const
{
  const double l = sqrt(z_shift * z_shift + R_shift * R_shift) + R_mirror_ref / 2;

  return l * sin(get_Rotation_HBD()) - get_RZ_Seg1_HBD() * cos(get_Rotation_HBD());
}

double RICH_Geometry::get_Z_Tip_HBD() const
{
  const double l = sqrt(z_shift * z_shift + R_shift * R_shift) + R_mirror_ref / 2;

  return l * cos(get_Rotation_HBD()) + get_RZ_Seg1_HBD() * sin(get_Rotation_HBD());
}

double RICH_Geometry::get_Rotation_HBD() const
{
  return atan(R_shift / z_shift);
}

int ePHENIXRICHConstruction::is_in_sector(G4VPhysicalVolume *volume) const
{
  if (sector_vec.find(volume) != sector_vec.end())
    return volume->GetCopyNo();
  else
    return -1;
}
