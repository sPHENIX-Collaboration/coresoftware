#include "PHG4ProjCrystalCalorimeterDetector.h"

#include "PHG4CrystalCalorimeterDetector.h"
#include "PHG4CrystalCalorimeterDisplayAction.h"

#include <Geant4/G4Cons.hh>
#include <Geant4/G4Element.hh>  // for G4Element
#include <Geant4/G4GenericTrap.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Trd.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>  // for vector

class G4VSolid;
class PHCompositeNode;
class PHG4CrystalCalorimeterSubsystem;

using namespace std;

//_______________________________________________________________________
PHG4ProjCrystalCalorimeterDetector::PHG4ProjCrystalCalorimeterDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, const std::string &dnam)
  : PHG4CrystalCalorimeterDetector(subsys, Node, dnam)
  ,
  //  _dx_front(50.19*mm),		//****************************************************************//
  //  _dy_front(50.19*mm),		//****************************************************************//
  //  _dx_back(59.3154545455*mm),		// PANDA eEMCAL Numbers: Crystals are 2.4cm * 2.4cm on front face //
  //  _dy_back(59.3154545455*mm),		//****************************************************************//
  //  _dz_crystal(90.000*mm),		//****************************************************************//
  _dx_front(41.44 * mm)
  , _dy_front(41.44 * mm)
  , _dx_back(48.97454545455 * mm)
  , _dy_back(48.97454545455 * mm)
  , _dz_crystal(90.000 * mm)
  , _crystallogicnameprefix("eEcalCrystal")
  , _4x4_construct_file("")
  , _overlapcheck_local(false)
{
}

//_______________________________________________________________________
PHG4ProjCrystalCalorimeterDetector::~PHG4ProjCrystalCalorimeterDetector()
{
}

//_______________________________________________________________________
int PHG4ProjCrystalCalorimeterDetector::IsInCrystalCalorimeter(G4VPhysicalVolume *volume) const
{
  // if hit is in absorber material
  //  bool isinabsorber = false;

  if (volume->GetName().find(_crystallogicnameprefix) != string::npos)
  {
    return 1;
  }
  else if (volume->GetName().find("arbon") != string::npos)
  {
    return -1;
  }

  return 0;
}

//_______________________________________________________________________
void PHG4ProjCrystalCalorimeterDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4ProjCrystalCalorimeterDetector: Begin Construction" << endl;
  }

  if (_mapping_tower_file.empty())
  {
    cout << "ERROR in PHG4ProjCrystalCalorimeterDetector: No tower mapping file specified. Abort detector construction." << endl;
    cout << "Please run SetTowerMappingFile( std::string filename ) first." << endl;
    exit(1);
  }

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  //G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4Material *Air = G4Material::GetMaterial("G4_Galactic");

  G4VSolid *ecal_envelope_cone = new G4Cons("eEcal_envelope_solid",
                                            _rMin1, _rMax1,
                                            _rMin2, _rMax2,
                                            _dZ / 2.,
                                            _sPhi, _dPhi);

  G4LogicalVolume *ecal_envelope_log = new G4LogicalVolume(ecal_envelope_cone, Air, G4String("eEcal_envelope"), 0, 0, 0);

  GetDisplayAction()->AddVolume(ecal_envelope_log, "Envelope");
  /* Define rotation attributes for envelope cone */
  G4RotationMatrix ecal_rotm;
  ecal_rotm.rotateX(_rot_in_x);
  ecal_rotm.rotateY(_rot_in_y);
  ecal_rotm.rotateZ(_rot_in_z);

  /* Place envelope cone in simulation */
  new G4PVPlacement(G4Transform3D(ecal_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z)),
                    ecal_envelope_log, "CrystalCalorimeter", logicWorld, 0, false, _overlapcheck_local);

  /* Construct crystal calorimeter within envelope */
  ConstructProjectiveCrystals(ecal_envelope_log);

  return;
}

void PHG4ProjCrystalCalorimeterDetector::GetCarbonFiberAdjustments(G4double &adjust_width, G4double &adjust_length)
{
  adjust_width = 0.1258525627 * mm;   //Because the crystals are slightly angled, the carbon fiber needs to be shortened
  adjust_length = 2.4824474402 * mm;  //      from the mother volume (to prevent clipping) by this amount.
}

void PHG4ProjCrystalCalorimeterDetector::GetCarbonFiberSpacing(G4double &CF_width, G4double &Air_CF, G4double &Air_Cry)
{
  //Parameters of the spacing given by PANDA document arXiv:0810.1216v1 Fig. 7.25
  CF_width = 0.18 * mm;  //Width of the carbon fiber which surrounds the crystal
  Air_CF = 0.24 * mm;    //Air gap between crystal and the carbon fiber
  Air_Cry = 0.60 * mm;   //Air gap between crystal and crystal
}

int PHG4ProjCrystalCalorimeterDetector::Fill4x4Unit(G4LogicalVolume *crystal_logic)
{
  //*************************************
  //**********Define Materials***********
  //*************************************

  //Crystal Material (Default is Lead Tungstate)
  G4Material *material_crystal = G4Material::GetMaterial(_materialCrystal.c_str());

  //Carbon Fiber
  G4double a = 12.01 * g / mole;
  G4Element *elC = new G4Element("Carbon", "C", 6., a);

  G4double density_carbon_fiber = 10 * 0.144 * g / cm3;
  G4Material *CarbonFiber = new G4Material("CarbonFiber", density_carbon_fiber, 1);
  CarbonFiber->AddElement(elC, 1);

  //Air
  //G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4Material *Air = G4Material::GetMaterial("G4_Galactic");

  //*************************************
  //**********Build First Crystal********
  //*************************************

  //Crystal Dimensions determined by the _dx_front, with various gaps and carbon fiber widths subtracted out

  //Parameters of the spacing given by PANDA document arXiv:0810.1216v1 Fig. 7.25
  G4double carbon_fiber_width, air_gap_carbon_fiber, air_gap_crystals;
  GetCarbonFiberSpacing(carbon_fiber_width, air_gap_carbon_fiber, air_gap_crystals);

  //Crystal Dimensions
  G4double dx_front_small = (_dx_front - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;  //Full width of the front crystal face
  G4double dy_front_small = (_dy_front - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;  //Full height of the front crystal face
  G4double dx_back_small = (_dx_back - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;    //Full width of the back crystal face
  G4double dy_back_small = (_dy_back - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;    //Full height of the back crystal face
  G4double dz = _dz_crystal;

  //Vertices of the primary, irregularly shaped crystal put into a vector
  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector(0, 0));
  vertices.push_back(G4TwoVector(0, dy_front_small));
  vertices.push_back(G4TwoVector(dx_front_small, dy_front_small));
  vertices.push_back(G4TwoVector(dx_front_small, 0));
  vertices.push_back(G4TwoVector(0, 0));
  vertices.push_back(G4TwoVector(0, dy_back_small));
  vertices.push_back(G4TwoVector(dx_back_small, dy_back_small));
  vertices.push_back(G4TwoVector(dx_back_small, 0));

  //Create the primary, irregularly shaped crystal
  G4VSolid *crystal_solid_small = new G4GenericTrap(G4String("eEcal_crystal"),
                                                    dz,
                                                    vertices);

  G4LogicalVolume *crystal_logic_small = new G4LogicalVolume(crystal_solid_small,
                                                             material_crystal,
                                                             "eEcal_crystal",
                                                             0, 0, 0);

  GetDisplayAction()->AddVolume(crystal_logic_small, "Crystal");

  //****************************************************
  //Build the solid/logical volume for the 2 x 2 crystal
  //****************************************************

  G4double TwoByTwo_dx1 = ((2.0 * dx_front_small) + (2.0 * air_gap_carbon_fiber) + air_gap_crystals) / 2.0;
  G4double TwoByTwo_dx2 = ((2.0 * dx_back_small) + (2.0 * air_gap_carbon_fiber) + air_gap_crystals) / 2.0;
  G4double TwoByTwo_dy1 = TwoByTwo_dx1;
  G4double TwoByTwo_dy2 = TwoByTwo_dx2;
  G4double TwoByTwo_dz = _dz_crystal;

  G4VSolid *Two_by_Two_solid = new G4Trd(G4String("Two_by_Two_solid"),
                                         TwoByTwo_dx1,  //Half length on the small face in x
                                         TwoByTwo_dx2,  //Half length on the large face in x
                                         TwoByTwo_dy1,  //Half length on the small face in y
                                         TwoByTwo_dy2,  //Half length on the large face in y
                                         TwoByTwo_dz);  //Half length in z

  G4LogicalVolume *Two_by_Two_logic = new G4LogicalVolume(Two_by_Two_solid,
                                                          Air,
                                                          "2_by_2_unit",
                                                          0, 0, 0);

  GetDisplayAction()->AddVolume(Two_by_Two_logic, "TwoByTwo");

  //*************************************************************
  //Read in mapping file for a single 2 x 2 block and 4 x 4 block
  //*************************************************************

  //The first four lines of the data file refer to the 2x2 block, and the last four lines refer to the mapping of the 4x4 block

  const string Crystal_Mapping_Small = _4x4_construct_file;  //Get the mapping file for the 4 x 4 block
  const int NumberOfIndices = 9;                             //Number of indices in mapping file for 4x4 block

  ifstream datafile_2;

  //Open the datafile, if it won't open return an error
  if (!datafile_2.is_open())
  {
    datafile_2.open(Crystal_Mapping_Small.c_str());
    if (!datafile_2)
    {
      cerr << endl
           << "*******************************************************************" << endl;
      cerr << "ERROR in 2 by 2 crystal mapping";
      cerr << "Failed to open " << Crystal_Mapping_Small << " --- Exiting program." << endl;
      cerr << "*******************************************************************" << endl
           << endl;
      exit(1);
    }
  }

  //Find the number of lines in the file, make and fill a NumberOfLines by NumberOfIndices matrix with contents of data file
  int NumberOfLines = 0;
  ifstream in(Crystal_Mapping_Small.c_str());
  std::string unused;
  while (std::getline(in, unused))
    ++NumberOfLines;

  G4int j_cry = NumberOfLines;
  G4int k_cry = NumberOfIndices;

  double TwoByTwo[j_cry][k_cry];

  G4int j = 0;
  G4int k = 0;

  while (j_cry > j)
  {
    while (k_cry > k)
    {
      datafile_2 >> TwoByTwo[j][k];
      k++;
    }
    j++;
    k = 0;
  }

  //**************************************************
  //Place the single crystal in the 2x2 volume 4 times
  //**************************************************

  G4int j_idx, k_idx;
  G4double x_cent, y_cent, z_cent, rot_x, rot_y, rot_z;
  G4int MappingIndex;

  j = 0;
  while (j_cry > j)
  {
    MappingIndex = TwoByTwo[j][8];
    if (MappingIndex == 1)
    {
      j_idx = TwoByTwo[j][0];
      k_idx = TwoByTwo[j][1];
      x_cent = TwoByTwo[j][2];
      y_cent = TwoByTwo[j][3];
      z_cent = TwoByTwo[j][4];
      //rot_x = TwoByTwo[j][5];
      //rot_y = TwoByTwo[j][6];
      rot_z = TwoByTwo[j][7];

      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(0 * rad);
      Rot->rotateY(0 * rad);
      Rot->rotateZ(rot_z * rad);

      ostringstream crystal_name;
      crystal_name.str("");
      crystal_name << _crystallogicnameprefix << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        crystal_logic_small,
                        crystal_name.str().c_str(),
                        Two_by_Two_logic,
                        0, 0, _overlapcheck_local);

      j_idx = k_idx = 0;
      x_cent = y_cent = z_cent = rot_x = rot_y = rot_z = 0.0;
    }
    j++;
  }

  //*************************************************************
  //Place the 2x2 volume in the 4x4 volume 4 times, with rotation
  //*************************************************************

  j = 0;
  while (j_cry > j)
  {
    MappingIndex = TwoByTwo[j][8];
    if (MappingIndex == 4)
    {
      j_idx = TwoByTwo[j][0];
      k_idx = TwoByTwo[j][1];
      x_cent = TwoByTwo[j][2];
      y_cent = TwoByTwo[j][3];
      z_cent = TwoByTwo[j][4];
      rot_x = TwoByTwo[j][5];
      rot_y = TwoByTwo[j][6];
      rot_z = TwoByTwo[j][7];

      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();
      Rot->rotateX(rot_x * rad);
      Rot->rotateY(rot_y * rad);
      Rot->rotateZ(0 * rad);

      ostringstream Two_by_Two_name;
      Two_by_Two_name.str("");
      Two_by_Two_name << "TwoByTwo"
                      << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        Two_by_Two_logic,
                        Two_by_Two_name.str().c_str(),
                        crystal_logic,
                        0, 0, _overlapcheck_local);

      j_idx = k_idx = 0;
      x_cent = y_cent = z_cent = rot_x = rot_y = 0.0;
    }

    j++;
  }

  //*****************************
  //Create the carbon fiber shell
  //*****************************

  //Create a hunk of carbon fiber the same size as the mother volume

  G4double dx1, dy1, dx2, dy2, dz_whole, carbon_fiber_adjust_width, carbon_fiber_adjust_length;

  GetCarbonFiberAdjustments(carbon_fiber_adjust_width, carbon_fiber_adjust_length);

  dx1 = _dx_front + carbon_fiber_adjust_width;
  dy1 = dx1;
  dx2 = _dx_back - carbon_fiber_adjust_width;
  dy2 = dx2;
  dz_whole = _dz_crystal - carbon_fiber_adjust_length;  //Carbon fiber shell should be slightly shorter than whole volume

  G4VSolid *Carbon_hunk_solid = new G4Trd(G4String("Carbon_hunk_solid"),
                                          dx1,        //Half length on the small face in x
                                          dx2,        //Half length on the large face in x
                                          dy1,        //Half length on the small face in y
                                          dy2,        //Half length on the large face in y
                                          dz_whole);  //Half length in z

  //Use subtraction solid to remove the crystal volumes from the carbon fiber hunk

  //----------------------------------------------------------------------------------------------------
  //First 2x2 crystal

  j = 0;
  G4int counter = 0;
  while (j_cry > counter)
  {
    MappingIndex = TwoByTwo[j][8];
    if (MappingIndex != 4) j++;
    counter++;
  }

  x_cent = TwoByTwo[j][2];
  y_cent = TwoByTwo[j][3];
  z_cent = TwoByTwo[j][4];
  rot_x = TwoByTwo[j][5];
  rot_y = TwoByTwo[j][6];
  //rot_z = TwoByTwo[j][7];

  G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

  G4RotationMatrix *Rot_1 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
  Rot_1->rotateX(rot_x * rad);
  Rot_1->rotateY(rot_y * rad);
  Rot_1->rotateZ(0 * rad);

  G4SubtractionSolid *Carbon_Shell_1 = new G4SubtractionSolid(G4String("Carbon_Shell_1"),
                                                              Carbon_hunk_solid,
                                                              Two_by_Two_solid,
                                                              Rot_1,
                                                              Crystal_Center);
  j++;

  //----------------------------------------------------------------------------------------------------
  //Second 2x2 crystal

  x_cent = TwoByTwo[j][2];
  y_cent = TwoByTwo[j][3];
  z_cent = TwoByTwo[j][4];
  rot_x = TwoByTwo[j][5];
  rot_y = TwoByTwo[j][6];
  //rot_z = TwoByTwo[j][7];

  Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

  G4RotationMatrix *Rot_2 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
  Rot_2->rotateX(rot_x * rad);
  Rot_2->rotateY(rot_y * rad);
  Rot_2->rotateZ(0 * rad);

  G4SubtractionSolid *Carbon_Shell_2 = new G4SubtractionSolid(G4String("Carbon_Shell_2"),
                                                              Carbon_Shell_1,
                                                              Two_by_Two_solid,
                                                              Rot_2,
                                                              Crystal_Center);
  j++;

  //----------------------------------------------------------------------------------------------------
  //Third 2x2 crystal

  x_cent = TwoByTwo[j][2];
  y_cent = TwoByTwo[j][3];
  z_cent = TwoByTwo[j][4];
  rot_x = TwoByTwo[j][5];
  rot_y = TwoByTwo[j][6];
  //rot_z = TwoByTwo[j][7];

  Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

  G4RotationMatrix *Rot_3 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
  Rot_3->rotateX(rot_x * rad);
  Rot_3->rotateY(rot_y * rad);
  Rot_3->rotateZ(0 * rad);

  G4SubtractionSolid *Carbon_Shell_3 = new G4SubtractionSolid(G4String("Carbon_Shell_3"),
                                                              Carbon_Shell_2,
                                                              Two_by_Two_solid,
                                                              Rot_3,
                                                              Crystal_Center);
  j++;

  //----------------------------------------------------------------------------------------------------
  //Final 2x2 crystal

  x_cent = TwoByTwo[j][2];
  y_cent = TwoByTwo[j][3];
  z_cent = TwoByTwo[j][4];
  rot_x = TwoByTwo[j][5];
  rot_y = TwoByTwo[j][6];
  //rot_z = TwoByTwo[j][7];

  Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

  G4RotationMatrix *Rot_4 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
  Rot_4->rotateX(rot_x * rad);
  Rot_4->rotateY(rot_y * rad);
  Rot_4->rotateZ(0 * rad);

  G4SubtractionSolid *Carbon_Shell_solid = new G4SubtractionSolid(G4String("Carbon_Shell_solid"),
                                                                  Carbon_Shell_3,
                                                                  Two_by_Two_solid,
                                                                  Rot_4,
                                                                  Crystal_Center);
  j++;

  //----------------------------------------------------------------------------------------------------

  //Create logical volume with the subtracted solid, made from carbon fiber material defined earlier

  G4LogicalVolume *Carbon_Shell_logic = new G4LogicalVolume(Carbon_Shell_solid,
                                                            CarbonFiber,
                                                            "Carbon_Fiber_logic",
                                                            0, 0, 0);

  GetDisplayAction()->AddVolume(Carbon_Shell_logic, "CarbonShell");

  //*************************************
  //Place the carbon fiber shell at (0,0)
  //*************************************

  G4ThreeVector Carbon_Center = G4ThreeVector(0 * mm, 0 * mm, 0 * mm);

  G4RotationMatrix *Rot_5 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
  Rot_5->rotateX(0 * rad);
  Rot_5->rotateY(0 * rad);
  Rot_5->rotateZ(0 * rad);

  new G4PVPlacement(Rot_5, Carbon_Center,
                    Carbon_Shell_logic,
                    "Carbon_Fiber_Shell",
                    crystal_logic,
                    0, 0, _overlapcheck_local);

  //***********************************
  //All done! Return to parent function
  //***********************************

  return 0;
}

int PHG4ProjCrystalCalorimeterDetector::FillSpecialUnit(G4LogicalVolume *crystal_logic, G4int ident)
{
  //*************************************
  //**********Define Materials***********
  //*************************************

  //Crystal Material (Default is Lead Tungstate)
  G4Material *material_crystal = G4Material::GetMaterial(_materialCrystal.c_str());

  //Carbon Fiber
  G4double a = 12.01 * g / mole;
  G4Element *elC = new G4Element("Carbon", "C", 6., a);

  G4double density_carbon_fiber = 10 * 0.144 * g / cm3;
  G4Material *CarbonFiber = new G4Material("CarbonFiber", density_carbon_fiber, 1);
  CarbonFiber->AddElement(elC, 1);

  //Air
  //G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4Material *Air = G4Material::GetMaterial("G4_Galactic");

  //*************************************
  //**********Build First Crystal********
  //*************************************

  //Crystal Dimensions determined by the _dx_front, with various gaps and carbon fiber widths subtracted out

  //Parameters of the spacing given by PANDA document arXiv:0810.1216v1 Fig. 7.25
  G4double carbon_fiber_width, air_gap_carbon_fiber, air_gap_crystals;
  GetCarbonFiberSpacing(carbon_fiber_width, air_gap_carbon_fiber, air_gap_crystals);

  //Crystal Dimensions
  G4double dx_front_small = (_dx_front - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;  //Full width of the front crystal face
  G4double dy_front_small = (_dy_front - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;  //Full height of the front crystal face
  G4double dx_back_small = (_dx_back - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;    //Full width of the back crystal face
  G4double dy_back_small = (_dy_back - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;    //Full height of the back crystal face
  G4double dz = _dz_crystal;

  //Vertices of the primary, irregularly shaped crystal put into a vector
  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector(0, 0));
  vertices.push_back(G4TwoVector(0, dy_front_small));
  vertices.push_back(G4TwoVector(dx_front_small, dy_front_small));
  vertices.push_back(G4TwoVector(dx_front_small, 0));
  vertices.push_back(G4TwoVector(0, 0));
  vertices.push_back(G4TwoVector(0, dy_back_small));
  vertices.push_back(G4TwoVector(dx_back_small, dy_back_small));
  vertices.push_back(G4TwoVector(dx_back_small, 0));

  //Create the primary, irregularly shaped crystal
  G4VSolid *crystal_solid_small = new G4GenericTrap(G4String("eEcal_crystal"),
                                                    dz,
                                                    vertices);

  G4LogicalVolume *crystal_logic_small = new G4LogicalVolume(crystal_solid_small,
                                                             material_crystal,
                                                             "eEcal_crystal",
                                                             0, 0, 0);

  GetDisplayAction()->AddVolume(crystal_logic_small, "Crystal");

  //****************************************************
  //Build the solid/logical volume for the 2 x 2 crystal
  //****************************************************

  G4double TwoByTwo_dx1 = ((2.0 * dx_front_small) + (2.0 * air_gap_carbon_fiber) + air_gap_crystals) / 2.0;
  G4double TwoByTwo_dx2 = ((2.0 * dx_back_small) + (2.0 * air_gap_carbon_fiber) + air_gap_crystals) / 2.0;
  G4double TwoByTwo_dy1 = TwoByTwo_dx1;
  G4double TwoByTwo_dy2 = TwoByTwo_dx2;
  G4double TwoByTwo_dz = _dz_crystal;

  G4VSolid *Two_by_Two_solid = new G4Trd(G4String("Two_by_Two_solid"),
                                         TwoByTwo_dx1,  //Half length on the small face in x
                                         TwoByTwo_dx2,  //Half length on the large face in x
                                         TwoByTwo_dy1,  //Half length on the small face in y
                                         TwoByTwo_dy2,  //Half length on the large face in y
                                         TwoByTwo_dz);  //Half length in z

  G4LogicalVolume *Two_by_Two_logic = new G4LogicalVolume(Two_by_Two_solid,
                                                          Air,
                                                          "2_by_2_unit",
                                                          0, 0, 0);

  GetDisplayAction()->AddVolume(Two_by_Two_logic, "TwoByTwo");

  //*************************************************************
  //Read in mapping file for a single 2 x 2 block and 4 x 4 block
  //*************************************************************

  //The first four lines of the data file refer to the 2x2 block, and the last four lines refer to the mapping of the 4x4 block

  const string Crystal_Mapping_Small = _4x4_construct_file;  //Get the mapping file for the 4 x 4 block
  const int NumberOfIndices = 9;                             //Number of indices in mapping file for 4x4 block

  ifstream datafile_2;

  //Open the datafile, if it won't open return an error
  if (!datafile_2.is_open())
  {
    datafile_2.open(Crystal_Mapping_Small.c_str());
    if (!datafile_2)
    {
      cerr << endl
           << "*******************************************************************" << endl;
      cerr << "ERROR in 2 by 2 crystal mapping";
      cerr << "Failed to open " << Crystal_Mapping_Small << " --- Exiting program." << endl;
      cerr << "*******************************************************************" << endl
           << endl;
      exit(1);
    }
  }

  //Find the number of lines in the file, make and fill a NumberOfLines by NumberOfIndices matrix with contents of data file
  int NumberOfLines = 0;
  ifstream in(Crystal_Mapping_Small.c_str());
  std::string unused;
  while (std::getline(in, unused))
    ++NumberOfLines;

  G4int j_cry = NumberOfLines;
  G4int k_cry = NumberOfIndices;

  double TwoByTwo[j_cry][k_cry];

  G4int j = 0;
  G4int k = 0;

  while (j_cry > j)
  {
    while (k_cry > k)
    {
      datafile_2 >> TwoByTwo[j][k];
      k++;
    }
    j++;
    k = 0;
  }

  //**************************************************
  //Place the single crystal in the 2x2 volume 4 times
  //**************************************************

  G4int j_idx, k_idx;
  G4double x_cent, y_cent, z_cent, rot_x, rot_y, rot_z;
  G4int MappingIndex;

  j = 0;
  while (j_cry > j)
  {
    MappingIndex = TwoByTwo[j][8];
    if (MappingIndex == 1)
    {
      j_idx = TwoByTwo[j][0];
      k_idx = TwoByTwo[j][1];
      x_cent = TwoByTwo[j][2];
      y_cent = TwoByTwo[j][3];
      z_cent = TwoByTwo[j][4];
      rot_x = TwoByTwo[j][5];
      rot_y = TwoByTwo[j][6];
      rot_z = TwoByTwo[j][7];

      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(0 * rad);
      Rot->rotateY(0 * rad);
      Rot->rotateZ(rot_z * rad);

      ostringstream crystal_name;
      crystal_name.str("");
      crystal_name << _crystallogicnameprefix << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        crystal_logic_small,
                        crystal_name.str().c_str(),
                        Two_by_Two_logic,
                        0, 0, _overlapcheck_local);

      j_idx = k_idx = 0;
      x_cent = y_cent = z_cent = rot_z = 0.0;
    }
    j++;
  }

  //****************************
  //Place the 2x2 Crystal Blocks
  //****************************

  j = 0;
  while (j_cry > j)
  {
    j_idx = TwoByTwo[j][0];
    k_idx = TwoByTwo[j][1];
    x_cent = TwoByTwo[j][2];
    y_cent = TwoByTwo[j][3];
    z_cent = TwoByTwo[j][4];
    rot_x = TwoByTwo[j][5];
    rot_y = TwoByTwo[j][6];
    rot_z = TwoByTwo[j][7];
    MappingIndex = TwoByTwo[j][8];

    if (ident == 12)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(0 * mm, 0 * mm, 0 * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();
      Rot->rotateX(0 * rad);
      Rot->rotateY(0 * rad);
      Rot->rotateZ(0 * rad);

      ostringstream Two_by_Two_name;
      Two_by_Two_name.str("");
      Two_by_Two_name << "TwoByTwo"
                      << "_j_" << 0 << "_k_" << 0;

      new G4PVPlacement(Rot, Crystal_Center,
                        Two_by_Two_logic,
                        Two_by_Two_name.str().c_str(),
                        crystal_logic,
                        0, 0, _overlapcheck_local);
    }
    else if (ident == 22)
    {
      if (MappingIndex == 22)
      {
        G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

        G4RotationMatrix *Rot = new G4RotationMatrix();
        Rot->rotateX(rot_x * rad);
        Rot->rotateY(0 * rad);
        Rot->rotateZ(0 * rad);

        ostringstream Two_by_Two_name;
        Two_by_Two_name.str("");
        Two_by_Two_name << "TwoByTwo"
                        << "_j_" << j_idx << "_k_" << k_idx;

        new G4PVPlacement(Rot, Crystal_Center,
                          Two_by_Two_logic,
                          Two_by_Two_name.str().c_str(),
                          crystal_logic,
                          0, 0, _overlapcheck_local);
      }
    }
    else if (ident == 32)
    {
      if (MappingIndex == 32)
      {
        G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

        G4RotationMatrix *Rot = new G4RotationMatrix();
        Rot->rotateX(rot_x * rad);
        Rot->rotateY(rot_y * rad);
        Rot->rotateZ(0 * rad);

        ostringstream Two_by_Two_name;
        Two_by_Two_name.str("");
        Two_by_Two_name << "TwoByTwo"
                        << "_j_" << j_idx << "_k_" << k_idx;

        new G4PVPlacement(Rot, Crystal_Center,
                          Two_by_Two_logic,
                          Two_by_Two_name.str().c_str(),
                          crystal_logic,
                          0, 0, _overlapcheck_local);
      }
    }
    else
    {
      cerr << endl
           << "Invalid Mapping Type: " << ident << endl;
      return -1;
    }

    j_idx = k_idx = 0;
    x_cent = y_cent = z_cent = rot_x = rot_y = rot_z = 0.0;
    j++;
  }

  //*******************
  //Create Carbon Fiber
  //*******************

  G4double dx1, dx2, dy1, dy2;
  GetCrystalSize(dx1, dy1, dx2, dy2, dz);

  if (ident == 12)
  {
    G4VSolid *Carbon_hunk_solid = new G4Trd(G4String("Carbon_hunk_solid"),
                                            (dx1 / 2.00 + 0.0209292929),  //Half length on the small face in x
                                            (dx2 / 2.00 - 0.0209292929),  //Half length on the large face in x
                                            (dy1 / 2.00 + 0.0209292929),  //Half length on the small face in y
                                            (dy2 / 2.0 - 0.0209292929),   //Half length on the large face in y
                                            (dz - 2 * mm));

    x_cent = 0 * mm;
    y_cent = 0 * mm;
    z_cent = 0 * mm;
    rot_x = 0 * mm;
    rot_y = 0 * mm;
    rot_z = 0 * mm;

    G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

    G4RotationMatrix *Rot_1 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
    Rot_1->rotateX(rot_x * rad);
    Rot_1->rotateY(rot_y * rad);
    Rot_1->rotateZ(rot_z * rad);

    G4SubtractionSolid *Carbon_Shell_solid = new G4SubtractionSolid(G4String("Carbon_Shell_solid"),
                                                                    Carbon_hunk_solid,
                                                                    Two_by_Two_solid,
                                                                    Rot_1,
                                                                    Crystal_Center);

    G4LogicalVolume *Carbon_Shell_logic = new G4LogicalVolume(Carbon_Shell_solid,
                                                              CarbonFiber,
                                                              "Carbon_Fiber_logic",
                                                              0, 0, 0);

    GetDisplayAction()->AddVolume(Carbon_Shell_logic, "CarbonShell");

    //*************************************
    //Place the carbon fiber shell at (0,0)
    //*************************************

    G4ThreeVector Carbon_Center = G4ThreeVector(0 * mm, 0 * mm, 0 * mm);

    G4RotationMatrix *Rot_5 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
    Rot_5->rotateX(0 * rad);
    Rot_5->rotateY(0 * rad);
    Rot_5->rotateZ(0 * rad);

    new G4PVPlacement(Rot_5, Carbon_Center,
                      Carbon_Shell_logic,
                      "Carbon_Fiber_Shell",
                      crystal_logic,
                      0, 0, _overlapcheck_local);
  }
  else if (ident == 22)
  {
    G4double carbon_fiber_adjust_width;   //Because the crystals are slightly angled, the carbon fiber needs to be shortened
    G4double carbon_fiber_adjust_length;  //	from the mother volume (to prevent clipping) by this amount.
    GetCarbonFiberAdjustments(carbon_fiber_adjust_width, carbon_fiber_adjust_length);
    G4double x_adjust = 0.0519558696 * mm;

    G4VSolid *Carbon_hunk_solid = new G4Trd(G4String("Carbon_hunk_solid"),
                                            (dx1 / 2.0 + (x_adjust)),            //Half length on the small face in x
                                            (dx2 / 2.0 - (x_adjust)),            //Half length on the large face in x
                                            (dy1 + carbon_fiber_adjust_width),   //Half length on the small face in y
                                            (dy2 - carbon_fiber_adjust_width),   //Half length on the large face in y
                                            (dz - carbon_fiber_adjust_length));  //Half length in z

    j = 0;
    G4int counter = 0;
    while (j_cry > counter)
    {
      MappingIndex = TwoByTwo[j][8];
      if (MappingIndex != 22) j++;
      counter++;
    }

    //First Hole

    x_cent = TwoByTwo[j][2];
    y_cent = TwoByTwo[j][3];
    z_cent = TwoByTwo[j][4];
    rot_x = TwoByTwo[j][5];
    rot_y = TwoByTwo[j][6];
    rot_z = TwoByTwo[j][7];

    G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

    G4RotationMatrix *Rot_1 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
    Rot_1->rotateX(rot_x * rad);
    Rot_1->rotateY(0 * rad);
    Rot_1->rotateZ(0 * rad);

    G4SubtractionSolid *Carbon_Shell_1 = new G4SubtractionSolid(G4String("Carbon_Shell_1"),
                                                                Carbon_hunk_solid,
                                                                Two_by_Two_solid,
                                                                Rot_1,
                                                                Crystal_Center);

    j++;

    //Second Hole

    x_cent = TwoByTwo[j][2];
    y_cent = TwoByTwo[j][3];
    z_cent = TwoByTwo[j][4];
    rot_x = TwoByTwo[j][5];
    //rot_y = TwoByTwo[j][6];
    //rot_z = TwoByTwo[j][7];

    Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

    G4RotationMatrix *Rot_2 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
    Rot_2->rotateX(rot_x * rad);
    Rot_2->rotateY(0 * rad);
    Rot_2->rotateZ(0 * rad);

    G4SubtractionSolid *Carbon_Shell_solid = new G4SubtractionSolid(G4String("Carbon_Shell_solid"),
                                                                    Carbon_Shell_1,
                                                                    Two_by_Two_solid,
                                                                    Rot_2,
                                                                    Crystal_Center);

    G4LogicalVolume *Carbon_Shell_logic = new G4LogicalVolume(Carbon_Shell_solid,
                                                              CarbonFiber,
                                                              "Carbon_Fiber_logic",
                                                              0, 0, 0);

    GetDisplayAction()->AddVolume(Carbon_Shell_logic, "CarbonShell");

    //*************************************
    //Place the carbon fiber shell at (0,0)
    //*************************************

    G4ThreeVector Carbon_Center = G4ThreeVector(0 * mm, 0 * mm, 0 * mm);

    G4RotationMatrix *Rot_5 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
    Rot_5->rotateX(0 * rad);
    Rot_5->rotateY(0 * rad);
    Rot_5->rotateZ(0 * rad);

    new G4PVPlacement(Rot_5, Carbon_Center,
                      Carbon_Shell_logic,
                      "Carbon_Fiber_Shell",
                      crystal_logic,
                      0, 0, _overlapcheck_local);
  }
  else if (ident == 32)
  {
    x_cent = 22.6036363636 + 0.18000 * 2.0;
    y_cent = x_cent;
    z_cent = 0.00;
    G4double rot_x = 0.020928529;
    G4double rot_y = -1.0 * rot_x;

    G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

    G4RotationMatrix *Rot_1 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
    Rot_1->rotateX(rot_x * rad);
    Rot_1->rotateY(rot_y * rad);
    Rot_1->rotateZ(0 * rad);

    G4VSolid *FourByFour_hunk_solid = new G4Trd(G4String("4x4_hunk_solid"),
                                                dx1,  //Half length on the small face in x
                                                dx2,  //Half length on the large face in x
                                                dy1,  //Half length on the small face in y
                                                dy2,  //Half length on the large face in y
                                                dz);  //Half length in z

    //Parameters of the spacing given by PANDA document arXiv:0810.1216v1 Fig. 7.25
    G4double carbon_fiber_width, air_gap_carbon_fiber, air_gap_crystals;
    GetCarbonFiberSpacing(carbon_fiber_width, air_gap_carbon_fiber, air_gap_crystals);

    //Crystal Dimensions
    G4double dx_front_small = (_dx_front - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;  //Full width of the front crystal face
    //G4double dy_front_small = ( _dy_front - (2.0 * carbon_fiber_width )  - (2.0 * air_gap_carbon_fiber) - air_gap_crystals ) / 2.0;		//Full height of the front crystal face
    G4double dx_back_small = (_dx_back - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;  //Full width of the back crystal face
    //G4double dy_back_small = (_dy_back - (2.0 * carbon_fiber_width )  - (2.0 * air_gap_carbon_fiber) - air_gap_crystals ) / 2.0;		//Full height of the back crystal face

    G4double TwoByTwo_dx1 = ((2.0 * dx_front_small) + (2.0 * air_gap_carbon_fiber) + air_gap_crystals) / 2.0;
    G4double TwoByTwo_dx2 = ((2.0 * dx_back_small) + (2.0 * air_gap_carbon_fiber) + air_gap_crystals) / 2.0;
    G4double TwoByTwo_dy1 = TwoByTwo_dx1;
    G4double TwoByTwo_dy2 = TwoByTwo_dx2;
    G4double TwoByTwo_dz = _dz_crystal;

    G4VSolid *Two_by_Two_solid = new G4Trd(G4String("Two_by_Two_solid"),
                                           TwoByTwo_dx1,           //Half length on the small face in x
                                           TwoByTwo_dx2,           //Half length on the large face in x
                                           TwoByTwo_dy1,           //Half length on the small face in y
                                           TwoByTwo_dy2,           //Half length on the large face in y
                                           TwoByTwo_dz + 2 * mm);  //Half length in z

    G4SubtractionSolid *Carbon_hunk_solid = new G4SubtractionSolid(G4String("Carbon_hunk_solid"),
                                                                   FourByFour_hunk_solid,
                                                                   Two_by_Two_solid,
                                                                   Rot_1,
                                                                   Crystal_Center);

    //----------------------------------------------------------------------------------------------------
    //First 2x2 crystal hole

    j = 0;
    G4int counter = 0;
    while (j_cry > counter)
    {
      MappingIndex = TwoByTwo[j][8];
      if (MappingIndex != 32) j++;
      counter++;
    }

    x_cent = TwoByTwo[j][2];
    y_cent = TwoByTwo[j][3];
    z_cent = TwoByTwo[j][4];
    rot_x = TwoByTwo[j][5];
    rot_y = TwoByTwo[j][6];
    rot_z = TwoByTwo[j][7];

    Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

    G4RotationMatrix *Rot_6 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
    Rot_6->rotateX(rot_x * rad);
    Rot_6->rotateY(rot_y * rad);
    Rot_6->rotateZ(0 * rad);

    G4SubtractionSolid *Carbon_Shell_1 = new G4SubtractionSolid(G4String("Carbon_Shell_1"),
                                                                Carbon_hunk_solid,
                                                                Two_by_Two_solid,
                                                                Rot_6,
                                                                Crystal_Center);
    j++;

    //----------------------------------------------------------------------------------------------------
    //Second 2x2 crystal hole

    x_cent = TwoByTwo[j][2];
    y_cent = TwoByTwo[j][3];
    z_cent = TwoByTwo[j][4];
    rot_x = TwoByTwo[j][5];
    rot_y = TwoByTwo[j][6];
    //rot_z = TwoByTwo[j][7];

    Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

    G4RotationMatrix *Rot_2 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
    Rot_2->rotateX(rot_x * rad);
    Rot_2->rotateY(rot_y * rad);
    Rot_2->rotateZ(0 * rad);

    G4SubtractionSolid *Carbon_Shell_2 = new G4SubtractionSolid(G4String("Carbon_Shell_2"),
                                                                Carbon_Shell_1,
                                                                Two_by_Two_solid,
                                                                Rot_2,
                                                                Crystal_Center);
    j++;

    //----------------------------------------------------------------------------------------------------
    //Third 2x2 crystal hole

    x_cent = TwoByTwo[j][2];
    y_cent = TwoByTwo[j][3];
    z_cent = TwoByTwo[j][4];
    rot_x = TwoByTwo[j][5];
    rot_y = TwoByTwo[j][6];
    //rot_z = TwoByTwo[j][7];

    Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

    G4RotationMatrix *Rot_3 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
    Rot_3->rotateX(rot_x * rad);
    Rot_3->rotateY(rot_y * rad);
    Rot_3->rotateZ(0 * rad);

    G4SubtractionSolid *Carbon_Shell_solid = new G4SubtractionSolid(G4String("Carbon_Shell_solid"),
                                                                    Carbon_Shell_2,
                                                                    Two_by_Two_solid,
                                                                    Rot_3,
                                                                    Crystal_Center);

    G4LogicalVolume *Carbon_Shell_logic = new G4LogicalVolume(Carbon_Shell_solid,
                                                              CarbonFiber,
                                                              "Carbon_Fiber_logic",
                                                              0, 0, 0);

    GetDisplayAction()->AddVolume(Carbon_Shell_logic, "CarbonShell");

    //*************************************
    //Place the carbon fiber shell at (0,0)
    //*************************************

    G4ThreeVector Carbon_Center = G4ThreeVector(0 * mm, 0 * mm, 0 * mm);

    G4RotationMatrix *Rot_5 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
    Rot_5->rotateX(0 * rad);
    Rot_5->rotateY(0 * rad);
    Rot_5->rotateZ(0 * rad);

    new G4PVPlacement(Rot_5, Carbon_Center,
                      Carbon_Shell_logic,
                      "Carbon_Fiber_Shell",
                      crystal_logic,
                      0, 0, _overlapcheck_local);
  }
  else
  {
    cerr << endl
         << "This is an error message which should never be shown. You may have disabled an important portion of this code upsream. Sorry :)" << endl;
    return -1;
  }

  return 0;
}

//_______________________________________________________________________
int PHG4ProjCrystalCalorimeterDetector::ConstructProjectiveCrystals(G4LogicalVolume *ecalenvelope)
{
  G4int NumberOfLines;                                  //Number of crystals to be created.
  const G4int NumberOfIndices = 9;                      //Different dimensions needed for crystal placement
  const string FileName = _mapping_tower_file.c_str();  //File in which crystal positions are stored

  G4int j_cry, k_cry;                               //Indices for matrix
  G4int j_idx, k_idx;                               //Indices of each crstals
  G4double x_cent, y_cent, z_cent, r_theta, r_phi;  //Coordinates of crystal in [x,y,z,theta,phi]

  G4double dx1;  //Half of the extent of the front face of the trapezoid in x
  G4double dx2;  //Half of the extent of the back face of the trapezoid in x
  G4double dy1;  //Half of the extent of the front face of the trapezoid in y
  G4double dy2;  //Half of the extent of the back face of the trapezoid in y
  G4double dz;   //Half of the extent of the crystal in z

  GetCrystalSize(dx1, dy1, dx2, dy2, dz);  //Fill crystal dimensions with function PHG4ProjCrystalCalorimeterDetector::CrystalDimensions

  //Create single crystal

  G4VSolid *crystal_solid = new G4Trd(G4String("eEcal_crystal"),
                                      dx1,  //Half length on the small face in x
                                      dx2,  //Half length on the large face in x
                                      dy1,  //Half length on the small face in y
                                      dy2,  //Half length on the large face in y
                                      dz);  //Half length in z

  //G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4Material *Air = G4Material::GetMaterial("G4_Galactic");

  G4LogicalVolume *crystal_logic = new G4LogicalVolume(crystal_solid,
                                                       Air,
                                                       "eEcal_crystal_unit",
                                                       0, 0, 0);

  //Create the single 2x2 unit (MappingIndex = 12)

  G4VSolid *twelve_solid = new G4Trd(G4String("12_solid"),
                                     (dx1 / 2.0),  //Half length on the small face in x
                                     (dx2 / 2.0),  //Half length on the large face in x
                                     (dy1 / 2.0),  //Half length on the small face in y
                                     (dy2 / 2.0),  //Half length on the large face in y
                                     dz);          //Half length in z

  G4LogicalVolume *twelve_logic = new G4LogicalVolume(twelve_solid,
                                                      Air,
                                                      "12_unit",
                                                      0, 0, 0);

  //Create the double 2x2 unit (i-shaped, MappingIndex = 22)

  G4VSolid *twentytwo_solid = new G4Trd(G4String("22_solid"),
                                        (dx1 / 2.0),  //Half length on the small face in x
                                        (dx2 / 2.0),  //Half length on the large face in x
                                        dy1,          //Half length on the small face in y
                                        dy2,          //Half length on the large face in y
                                        dz);          //Half length in z

  G4LogicalVolume *twentytwo_logic = new G4LogicalVolume(twentytwo_solid,
                                                         Air,
                                                         "22_unit",
                                                         0, 0, 0);

  //Create the triple 2x2 unit (L-shaped, MappingIndex = 32)

  x_cent = 22.6036363636 + 0.18000 * 2.0;
  y_cent = x_cent;
  z_cent = 0.00;
  G4double rot_x = 0.020928529;
  G4double rot_y = -1.0 * rot_x;
  G4double rot_z = 0.00;

  G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

  G4RotationMatrix *Rot_1 = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
  Rot_1->rotateX(rot_x * rad);
  Rot_1->rotateY(rot_y * rad);
  Rot_1->rotateZ(0 * rad);

  G4VSolid *FourByFour_hunk_solid = new G4Trd(G4String("4x4_hunk_solid"),
                                              dx1,  //Half length on the small face in x
                                              dx2,  //Half length on the large face in x
                                              dy1,  //Half length on the small face in y
                                              dy2,  //Half length on the large face in y
                                              dz);  //Half length in z

  //Parameters of the spacing given by PANDA document arXiv:0810.1216v1 Fig. 7.25
  G4double carbon_fiber_width = 0.18 * mm;    //Width of the carbon fiber which surrounds the crystal
  G4double air_gap_carbon_fiber = 0.24 * mm;  //Air gap between crystal and the carbon fiber
  G4double air_gap_crystals = 0.60 * mm;      //Air gap between crystal and crystal

  //Crystal Dimensions
  G4double dx_front_small = (_dx_front - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;  //Full width of the front crystal face
                                                                                                                               //	G4double dy_front_small = ( _dy_front - (2.0 * carbon_fiber_width )  - (2.0 * air_gap_carbon_fiber) - air_gap_crystals ) / 2.0;		//Full height of the front crystal face
  G4double dx_back_small = (_dx_back - (2.0 * carbon_fiber_width) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals) / 2.0;    //Full width of the back crystal face
                                                                                                                               //	G4double dy_back_small = (_dy_back - (2.0 * carbon_fiber_width )  - (2.0 * air_gap_carbon_fiber) - air_gap_crystals ) / 2.0;		//Full height of the back crystal face

  G4double TwoByTwo_dx1 = ((2.0 * dx_front_small) + (2.0 * air_gap_carbon_fiber) + air_gap_crystals) / 2.0;
  G4double TwoByTwo_dx2 = ((2.0 * dx_back_small) + (2.0 * air_gap_carbon_fiber) + air_gap_crystals) / 2.0;
  G4double TwoByTwo_dy1 = TwoByTwo_dx1;
  G4double TwoByTwo_dy2 = TwoByTwo_dx2;
  G4double TwoByTwo_dz = _dz_crystal;

  G4VSolid *Two_by_Two_solid = new G4Trd(G4String("Two_by_Two_solid"),
                                         TwoByTwo_dx1,           //Half length on the small face in x
                                         TwoByTwo_dx2,           //Half length on the large face in x
                                         TwoByTwo_dy1,           //Half length on the small face in y
                                         TwoByTwo_dy2,           //Half length on the large face in y
                                         TwoByTwo_dz + 5 * mm);  //Half length in z

  G4SubtractionSolid *thirtytwo_solid = new G4SubtractionSolid(G4String("32_solid"),
                                                               FourByFour_hunk_solid,
                                                               Two_by_Two_solid,
                                                               Rot_1,
                                                               Crystal_Center);

  G4LogicalVolume *thirtytwo_logic = new G4LogicalVolume(thirtytwo_solid,
                                                         Air,
                                                         "32_unit",
                                                         0, 0, 0);

  GetDisplayAction()->AddVolume(crystal_logic, "Invisible");
  GetDisplayAction()->AddVolume(twelve_logic, "Invisible");
  GetDisplayAction()->AddVolume(twentytwo_logic, "Invisible");
  GetDisplayAction()->AddVolume(thirtytwo_logic, "Invisible");

  Fill4x4Unit(crystal_logic);
  FillSpecialUnit(twelve_logic, 12);
  FillSpecialUnit(twentytwo_logic, 22);
  FillSpecialUnit(thirtytwo_logic, 32);

  ostringstream name;

  ifstream datafile;

  if (!datafile.is_open())
  {
    datafile.open(FileName.c_str());
    if (!datafile)
    {
      cerr << endl
           << "*******************************************************************" << endl;
      cerr << "ERROR: Failed to open " << FileName << " --- Exiting program." << endl;
      cerr << "*******************************************************************" << endl
           << endl;
      exit(1);
    }
  }

  //Determine the number of crystals to be created
  NumberOfLines = 0;
  ifstream in(_mapping_tower_file.c_str());
  std::string unused;
  while (std::getline(in, unused))
    ++NumberOfLines;

  j_cry = NumberOfLines;    // = Number of Crystals
  k_cry = NumberOfIndices;  // = j, k, x, y, z, alpha, beta.

  double Crystals[j_cry][k_cry];

  G4int j = 0;
  G4int k = 0;

  //Fill matrix with the data from the mapping file
  while (j_cry > j)
  {
    while (k_cry > k)
    {
      datafile >> Crystals[j][k];
      k++;
    }
    j++;
    k = 0;
  }

  // Find the maximum k index
  G4int k_max = 0;
  G4int MappingIndex;
  j = 0;
  while (j_cry > j)
  {
    MappingIndex = Crystals[j][8];
    if (Crystals[j][1] > k_max) k_max = Crystals[j][1];
    j++;
  }

  //Second Quadrant
  j = 0;
  while (j_cry > j)
  {
    MappingIndex = Crystals[j][8];

    j_idx = Crystals[j][0];
    k_idx = Crystals[j][1];
    x_cent = Crystals[j][2] - _place_in_x;
    y_cent = Crystals[j][3] - _place_in_y;
    z_cent = Crystals[j][4] - _place_in_z;  //Coordinate system refers to mother volume, have to subtract out its position in the actual xyz-space
    r_theta = Crystals[j][5];               //Rotation in Horizontal
    r_phi = Crystals[j][6];                 //Rotation in Vertical
    rot_z = Crystals[j][7];

    if (MappingIndex == 16)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(0 * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        crystal_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 32)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      if (j_idx == (k_max - 1) / 2 && k_idx == (k_max + 1) / 2) rot_z = rot_z / 2.0;

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        thirtytwo_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 22)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        twentytwo_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 12)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        twelve_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else
    {
    }

    j++;
  }

  //First Quadrant
  j = 0;
  while (j_cry > j)
  {
    MappingIndex = Crystals[j][8];

    j_idx = k_max - Crystals[j][0];
    k_idx = Crystals[j][1];
    x_cent = -1.0 * (Crystals[j][2] - _place_in_x);
    y_cent = Crystals[j][3] - _place_in_y;
    z_cent = Crystals[j][4] - _place_in_z;
    r_theta = -1.0 * Crystals[j][5];
    r_phi = Crystals[j][6];
    rot_z = Crystals[j][7] / (-4.0);

    if (MappingIndex == 16)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(0 * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        crystal_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 32)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      if (j_idx == (k_max + 1) / 2 && k_idx == (k_max + 1) / 2) rot_z = rot_z * 3.0;

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        thirtytwo_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 22)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      rot_z = Crystals[j][7] * (-1.0);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        twentytwo_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 12)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      rot_z = rot_z * -4.0000;

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        twelve_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else
    {
    }

    j++;
  }

  //Fourth Quadrant
  j = 0;
  while (j_cry > j)
  {
    MappingIndex = Crystals[j][8];
    j_idx = k_max - Crystals[j][0];
    k_idx = k_max - Crystals[j][1];
    x_cent = -1.0 * (Crystals[j][2] - _place_in_x);
    y_cent = -1.0 * (Crystals[j][3] - _place_in_y);
    z_cent = Crystals[j][4] - _place_in_z;
    r_theta = -1.0 * Crystals[j][5];
    r_phi = -1.0 * Crystals[j][6];
    rot_z = Crystals[j][7] / 2.0;

    if (MappingIndex == 16)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(0 * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        crystal_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 32)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      if (j_idx == (k_max + 1) / 2 && k_idx == (k_max - 1) / 2) rot_z = rot_z * -2.0;

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        thirtytwo_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 22)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      rot_z = Crystals[j][7] * (-1.0);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        twentytwo_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 12)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        twelve_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else
    {
    }

    j++;
  }

  //Third Quadrant
  j = 0;
  while (j_cry > j)
  {
    MappingIndex = Crystals[j][8];

    j_idx = Crystals[j][0];
    k_idx = k_max - Crystals[j][1];
    x_cent = Crystals[j][2] - _place_in_x;
    y_cent = -1.0 * (Crystals[j][3] - _place_in_y);
    z_cent = Crystals[j][4] - _place_in_z;
    r_theta = Crystals[j][5];
    r_phi = -1.0 * Crystals[j][6];
    rot_z = Crystals[j][7] / 4.0;

    if (MappingIndex == 16)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(0 * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        crystal_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 32)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      if (j_idx == (k_max - 1) / 2 && k_idx == (k_max - 1) / 2) rot_z = rot_z * 3.0;

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        thirtytwo_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 22)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      rot_z = Crystals[j][7] * (1.0);

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        twentytwo_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else if (MappingIndex == 12)
    {
      G4ThreeVector Crystal_Center = G4ThreeVector(x_cent * mm, y_cent * mm, z_cent * mm);

      rot_z = rot_z * -4.0000;

      G4RotationMatrix *Rot = new G4RotationMatrix();  //rotation matrix for the placement of each crystal
      Rot->rotateX(r_phi * rad);
      Rot->rotateY(r_theta * rad);
      Rot->rotateZ(rot_z * rad);

      name.str("");
      name << "FourByFour"
           << "_j_" << j_idx << "_k_" << k_idx;

      new G4PVPlacement(Rot, Crystal_Center,
                        twelve_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, _overlapcheck_local);
    }
    else
    {
    }

    j++;
  }

  return 0;
}
