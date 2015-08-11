#include "PHG4CrystalCalorimeterDetector.h"
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
#include <Geant4/G4GenericTrap.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Trd.hh>
#include <Geant4/G4Tubs.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <cmath>
#include <sstream>

#include <iostream>
#include <fstream>
#include <cstdlib>


using namespace std;

//static double no_overlap = 0.00015 * cm; // added safety margin against overlaps by using same boundary between volumes
static bool overlapcheck_local = true;

//_______________________________________________________________________
PHG4CrystalCalorimeterDetector::PHG4CrystalCalorimeterDetector( PHCompositeNode *Node, const std::string &dnam ):
  PHG4Detector(Node, dnam),
  _place_in_x(0.0*mm),
  _place_in_y(0.0*mm),
  _place_in_z(-1080.0*mm),
  _rot_in_x(0.0),
  _rot_in_y(M_PI),
  _rot_in_z(0.0),
  _rMin1(22*mm),
  _rMax1(656*mm),
  _rMin2(26*mm),
  _rMax2(775*mm),
  _dZ(180*mm),
  _sPhi(0),
  _dPhi(2*M_PI),
  _dx_front(41.44*mm),
  _dy_front(41.44*mm),
  _dx_back(48.97454545455*mm),
  _dy_back(48.97454545455*mm),
  _dz_crystal(90.000*mm),
  _materialCrystal( "G4_PbWO4" ),
  _active(1),
  _crystallogicnameprefix("eEcalCrystal"),
  _superdetector("NONE"),
  _inputFile( "" )
{

}


//_______________________________________________________________________
PHG4CrystalCalorimeterDetector::~PHG4CrystalCalorimeterDetector()
{}


//_______________________________________________________________________
int
PHG4CrystalCalorimeterDetector::IsInCrystalCalorimeter(G4VPhysicalVolume * volume) const
{
  // if hit is in absorber material
  //  bool isinabsorber = false;

  if (volume->GetName().find(_crystallogicnameprefix) != string::npos)
    {
      return 1;
    }
  else if ( volume->GetName().find("arbon") != string::npos )
    {
      return -1;
    }

  return 0;

}


//_______________________________________________________________________
void
PHG4CrystalCalorimeterDetector::Construct( G4LogicalVolume* logicWorld )
{
  if ( verbosity > 0 )
    {
      cout << "PHG4CrystalCalorimeterDetector: Begin Construction" << endl;
    }


  if ( _inputFile.empty() )
    {
      cout << "ERROR in PHG4CrystalCalorimeterDetector: No tower mapping file specified. Abort detector construction." << endl;
      cout << "Please run SetTowerMappingFile( std::string filename ) first." << endl;
      exit(1);
    }

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  //G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4Material* Air = G4Material::GetMaterial("G4_Galactic");

  G4VSolid* ecal_envelope_cone = new G4Cons("eEcal_envelope_solid",
					    _rMin1, _rMax1,
					    _rMin2, _rMax2,
					    _dZ/2.,
					    _sPhi, _dPhi );

  G4LogicalVolume* ecal_envelope_log =  new G4LogicalVolume(ecal_envelope_cone, Air, G4String("eEcal_envelope"), 0, 0, 0);

  /* Define visualization attributes for envelope cone */
  G4VisAttributes* ecalVisAtt = new G4VisAttributes();
	  ecalVisAtt->SetVisibility(false);
	  ecalVisAtt->SetForceSolid(false);
	  ecalVisAtt->SetColour(G4Colour::Magenta());
  ecal_envelope_log->SetVisAttributes(ecalVisAtt);

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix ecal_rotm;
  ecal_rotm.rotateX(_rot_in_x);
  ecal_rotm.rotateY(_rot_in_y);
  ecal_rotm.rotateZ(_rot_in_z);

  /* Place envelope cone in simulation */
  new G4PVPlacement( G4Transform3D(ecal_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z) ),
		     ecal_envelope_log, "CrystalCalorimeter", logicWorld, 0, false, overlapcheck_local);

  /* Construct crystal calorimeter within envelope */
  ConstructCrystals(ecal_envelope_log);

  return;
}

void
PHG4CrystalCalorimeterDetector::CrystalDimensions(G4double& dx_front, G4double& dy_front, G4double& dx_back, G4double& dy_back, G4double& dz)
{
	dx_front = _dx_front;
	dy_front = _dy_front;
	dx_back = _dx_back;
	dy_back = _dy_back;
	dz = _dz_crystal;

}


void
PHG4CrystalCalorimeterDetector::CarbonFiberAdjustments(G4double& adjust_width, G4double& adjust_length)
{
        adjust_width = 0.1258525627*mm;   //Because the crystals are slightly angled, the carbon fiber needs to be shortened 
        adjust_length = 2.4824474402*mm;  //      from the mother volume (to prevent clipping) by this amount.
}


void
PHG4CrystalCalorimeterDetector::CarbonFiberSpacing(G4double& CF_width, G4double& Air_CF, G4double& Air_Cry)
{
        //Parameters of the spacing given by PANDA document arXiv:0810.1216v1 Fig. 7.25
        CF_width = 0.18*mm;		//Width of the carbon fiber which surrounds the crystal
        Air_CF = 0.24*mm;		//Air gap between crystal and the carbon fiber
        Air_Cry = 0.60*mm;		//Air gap between crystal and crystal
}


int
PHG4CrystalCalorimeterDetector::ConstructCrystals(G4LogicalVolume* ecalenvelope)
{
  /* Simple model: Fill calorimeter envelope with rectangular crystals (front face = back face)
     which are arranged in chessboard pattern in x-y and oriented parallel to z azis */

  G4Material* Vacuum = G4Material::GetMaterial("G4_Galactic");
  //G4Material* material_crystal = G4Material::GetMaterial(_materialCrystal.c_str());

  G4double carbonThickness, airGap, airGap_Crystal;
  CarbonFiberSpacing( carbonThickness, airGap, airGap_Crystal);

  G4double crystal_dx = 20*mm + carbonThickness;
  G4double crystal_dy = 20*mm + carbonThickness;
  G4double crystal_dz = _dz_crystal;

  G4VSolid* crystal_solid = new G4Box( G4String("crystal_unit_solid"),
				       crystal_dx / 2.0,
				       crystal_dy / 2.0,
				       crystal_dz );

  G4LogicalVolume *crystal_logic = new G4LogicalVolume( crystal_solid,
							Vacuum,
							"crystal_unit_logical",
							0, 0, 0);

  FillTower(crystal_logic);

  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(false);
  visattchk->SetForceSolid(true);
  visattchk->SetColour(G4Colour::Yellow());
  crystal_logic->SetVisAttributes(visattchk);

  /* Place crystal units / tower */
  ifstream istream_mapping;

  /* Open the datafile, if it won't open return an error */
  if (!istream_mapping.is_open())
    {
      istream_mapping.open( _inputFile.c_str() );
      if(!istream_mapping)
	{
	  cerr << "ERROR in PHG4CrystalCalorimeterDetector: Failed to open mapping file " << _inputFile << endl;
	  exit(1);
	}
    }

  string line_mapping;

  while ( getline( istream_mapping, line_mapping ) )
    {
      unsigned idx_j, idx_k, idx_l;
      double pos_x, pos_y, pos_z;
      double dummy;

      istringstream iss(line_mapping);

      /* Skip lines starting with / including a '#' */
      if ( line_mapping.find("#") != string::npos )
	{
	  if ( verbosity > 0 )
	    {
	      cout << "PHG4CrystalCalorimeterDetector: SKIPPING line in mapping file: " << line_mapping << endl;
	    }
	  continue;
	}

      /* read string- break if error */
      if ( !( iss >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> dummy >> dummy >> dummy >> dummy ) )
	{
	  cerr << "ERROR in PHG4CrystalCalorimeterDetector: Failed to read line in mapping file " << _inputFile << endl;
	  exit(1);
	}

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      ostringstream towername;
      towername.str("");
      towername << _crystallogicnameprefix << "_j_" << idx_j << "_k_" << idx_k;

      /* Place tower */
      if ( verbosity > 0 )
	{
	  cout << "PHG4CrystalCalorimeterDetector: Place tower " << towername.str() << endl;
	}

      new G4PVPlacement( 0, G4ThreeVector(pos_x*cm , pos_y*cm, pos_z*cm),
			 crystal_logic,
			 towername.str().c_str(),
			 ecalenvelope,
			 0, 0, overlapcheck);
    }

  return 0;

}


int
PHG4CrystalCalorimeterDetector::FillTower(G4LogicalVolume *crystal_logic)
{

	//*************************************
	//**********Define Materials***********
	//*************************************
	
	//Crystal Material (Default is Lead Tungstate)
	G4Material* material_crystal = G4Material::GetMaterial(_materialCrystal.c_str());

	//Carbon Fiber
	G4double a = 12.01*g/mole;
	G4Element* elC = new G4Element("Carbon", "C", 6., a);
	
	G4double density_carbon_fiber = 0.144*g/cm3;
	G4Material* CarbonFiber = new G4Material("CarbonFiber", density_carbon_fiber, 1);
		CarbonFiber->AddElement(elC, 1);

	//Vacuum
	//G4Material* Vacuum = G4Material::GetMaterial("G4_AIR");

	//*************************************
	//**********Define Constants***********
	//*************************************
		
	G4double carbonThickness, airGap, airGap_Crystal;
	CarbonFiberSpacing( carbonThickness, airGap, airGap_Crystal);
	
	const G4double crystal_dx = 20*mm + carbonThickness; 
	const G4double crystal_dy = 20*mm + carbonThickness; 
	const G4double crystal_dz =  _dz_crystal;

	//*************************************
	//**********Build First Crystal********
	//*************************************

	ostringstream crystal_solid_name;
	crystal_solid_name.str("");
	crystal_solid_name << _crystallogicnameprefix << "_solid"; 

	G4double lead_dx = (( crystal_dx / 2.0 ) - carbonThickness);
	G4double lead_dy = (( crystal_dy / 2.0 ) - carbonThickness);
	
	G4VSolid* lead_solid = new G4Box( crystal_solid_name.str().c_str(),
		lead_dx,
		lead_dy,
		crystal_dz );

	ostringstream crystal_logic_name;
	crystal_logic_name.str("");
	crystal_solid_name << _crystallogicnameprefix << "_logic"; 
	
	G4LogicalVolume *lead_logic = new G4LogicalVolume( lead_solid,
		material_crystal,
		crystal_logic_name.str().c_str(),
		0, 0, 0);
	
	G4VisAttributes *visattchk = new G4VisAttributes();
		visattchk->SetVisibility(true);
		visattchk->SetForceSolid(true);
		visattchk->SetColour(G4Colour::Cyan());
	lead_logic->SetVisAttributes(visattchk);
	
	//********************************
	//Place the eECAL Crystal at (0,0)
	//********************************

	G4ThreeVector Crystal_Center = G4ThreeVector(0.00*mm, 0.00*mm, 0.00*mm);

	G4RotationMatrix *Rot = new G4RotationMatrix(); //rotation matrix for the placement of each crystal
		Rot->rotateX(0*rad);
		Rot->rotateY(0*rad);
		Rot->rotateZ(0*rad);

	ostringstream crystal_name;
	crystal_name.str("");
	crystal_name << _crystallogicnameprefix << "_default";
	
	new G4PVPlacement( Rot, Crystal_Center,
		lead_logic,
		crystal_name.str().c_str(),
		crystal_logic,
		0, 0, overlapcheck);
	
	//*****************************
	//Create the carbon fiber shell
	//*****************************

	G4VSolid* Carbon_hunk_solid = new G4Box( G4String("Carbon_hunk_solid"),
		crystal_dx / 2.0,
		crystal_dy / 2.0,
		( (crystal_dz) - 1*mm) );
	
	G4SubtractionSolid* Carbon_Shell = new G4SubtractionSolid(G4String("Carbon_Shell_solid"),
		Carbon_hunk_solid,
		lead_solid,
		Rot,
		Crystal_Center);

	G4LogicalVolume *Carbon_Shell_logic = new G4LogicalVolume( Carbon_Shell,
		CarbonFiber,
		G4String("Carbon_Shell"),
		0, 0, 0);

	G4VisAttributes *visattchk2 = new G4VisAttributes();
		visattchk2->SetVisibility(true);
		visattchk2->SetForceSolid(true);
		visattchk2->SetColour(G4Colour::Black());
	Carbon_Shell_logic->SetVisAttributes(visattchk2);

	//*************************************
	//Place the carbon fiber shell at (0,0)
	//*************************************

	new G4PVPlacement( Rot, Crystal_Center,
		Carbon_Shell_logic,
		G4String("Carbon_Shell"),
		crystal_logic,
		0, 0, overlapcheck);
	
	return 0;
}
