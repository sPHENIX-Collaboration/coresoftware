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

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <cmath>
#include <sstream>

#include <iostream>
#include <fstream>
#include <cstdlib>


using namespace std;


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
  _crystal_dx(20*mm),
  _crystal_dy(20*mm),
  _crystal_dz(180.0*mm),
  _materialCrystal( "G4_PbWO4" ),
  _active(1),
  _towerlogicnameprefix("CrystalCalorimeterTower"),
  _superdetector("NONE"),
  _mapping_tower_file("")
{

}


//_______________________________________________________________________
PHG4CrystalCalorimeterDetector::~PHG4CrystalCalorimeterDetector()
{}


//_______________________________________________________________________
int
PHG4CrystalCalorimeterDetector::IsInCrystalCalorimeter(G4VPhysicalVolume * volume) const
{

  if (volume->GetName().find(_towerlogicnameprefix) != string::npos)
    {
      if (volume->GetName().find("crystal") != string::npos)
	{
	  return 1;
	}
      /* only record energy in actual absorber- drop energy lost in air gaps inside envelope */
      else if (volume->GetName().find("absorber") != string::npos)
	{
	  return -1;
	}
      else if (volume->GetName().find("envelope") != string::npos)
	{
	  return 0;
	}
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


  if ( _mapping_tower_file.empty() )
    {
      cout << "ERROR in PHG4CrystalCalorimeterDetector: No tower mapping file specified. Abort detector construction." << endl;
      cout << "Please run SetTowerMappingFile( std::string filename ) first." << endl;
      exit(1);
    }

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* eemc_envelope_solid = new G4Cons("eemc_envelope_solid",
					    _rMin1, _rMax1,
					    _rMin2, _rMax2,
					    _dZ/2.,
					    _sPhi, _dPhi );

  G4LogicalVolume* eemc_envelope_log =  new G4LogicalVolume(eemc_envelope_solid, Air, G4String("eemc_envelope"), 0, 0, 0);

  /* Define visualization attributes for envelope cone */
  G4VisAttributes* eemcVisAtt = new G4VisAttributes();
  eemcVisAtt->SetVisibility(false);
  eemcVisAtt->SetForceSolid(false);
  eemcVisAtt->SetColour(G4Colour::Magenta());
  eemc_envelope_log->SetVisAttributes(eemcVisAtt);

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix eemc_rotm;
  eemc_rotm.rotateX(_rot_in_x);
  eemc_rotm.rotateY(_rot_in_y);
  eemc_rotm.rotateZ(_rot_in_z);

  /* Place envelope cone in simulation */
  ostringstream name_envelope;
  name_envelope.str("");
  name_envelope << _towerlogicnameprefix << "_envelope" << endl;

  new G4PVPlacement( G4Transform3D(eemc_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z) ),
		     eemc_envelope_log, name_envelope.str().c_str(), logicWorld, 0, false, overlapcheck);

  /* Construct single calorimeter tower */
  G4LogicalVolume* singletower = ConstructTower();

  /* Place calorimeter tower within envelope */
  PlaceTower( eemc_envelope_log , singletower );

  return;
}


//_______________________________________________________________________
G4LogicalVolume*
PHG4CrystalCalorimeterDetector::ConstructTower()
{
  if ( verbosity > 0 )
    {
      cout << "PHG4CrystalCalorimeterDetector: Build logical volume for single tower..." << endl;
    }

  G4double carbon_thickness = 0.009 * cm;
  G4double airgap_crystal_carbon = 0.012 * cm;

  /* dimesnions of full tower */
  G4double tower_dx = _crystal_dx + 2 * ( carbon_thickness + airgap_crystal_carbon );
  G4double tower_dy = _crystal_dy + 2 * ( carbon_thickness + airgap_crystal_carbon );
  G4double tower_dz = _crystal_dz;

  /* create logical volume for single tower */
  G4Material* material_air = G4Material::GetMaterial( "G4_AIR" );

  G4VSolid* single_tower_solid = new G4Box( G4String("single_tower_solid"),
					    tower_dx / 2.0,
					    tower_dy / 2.0,
					    tower_dz / 2.0 );

  G4LogicalVolume *single_tower_logic = new G4LogicalVolume( single_tower_solid,
							     material_air,
							     "single_tower_logic",
							     0, 0, 0);

  /* create geometry volume for crystal inside single_tower */
  G4VSolid* solid_crystal = new G4Box( G4String("single_crystal_solid"),
				       _crystal_dx / 2.0,
				       _crystal_dy / 2.0,
				       _crystal_dz / 2.0 );

  /* create geometry volume for frame (carbon fober shell) inside single_tower */
  G4VSolid* Carbon_hunk_solid = new G4Box( G4String("Carbon_hunk_solid"),
					   tower_dx / 2.0,
					   tower_dy / 2.0,
					   ( (tower_dz / 2.0) - 1*mm) );

  G4double lead_dx = tower_dx - 2.0 * carbon_thickness;
  G4double lead_dy = tower_dy - 2.0 * carbon_thickness;

  G4VSolid* lead_solid = new G4Box( G4String("lead_solid"),
				    lead_dx / 2.0,
				    lead_dy / 2.0,
				    tower_dz / 2.0 );

  G4SubtractionSolid* Carbon_shell_solid = new G4SubtractionSolid(G4String("Carbon_Shell_solid"),
								  Carbon_hunk_solid,
								  lead_solid,
								  0,
								  G4ThreeVector(0.00*mm, 0.00*mm, 0.00*mm) );


  /* create logical volumes for crystal inside single_tower */
  G4Material* material_crystal = G4Material::GetMaterial( _materialCrystal.c_str() );

  G4LogicalVolume *logic_crystal = new G4LogicalVolume( solid_crystal,
							material_crystal,
							"single_crystal_logic",
							0, 0, 0);

  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(true);
  visattchk->SetColour(G4Colour::Cyan());
  logic_crystal->SetVisAttributes(visattchk);

  /* create logical volumes for structural frame */
  //Carbon Fiber
  G4double a = 12.01*g/mole;
  G4Element* elC = new G4Element("Carbon", "C", 6., a);
  G4double density_carbon_fiber = 0.144*g/cm3;
  G4Material* CarbonFiber = new G4Material("CarbonFiber", density_carbon_fiber, 1);
  CarbonFiber->AddElement(elC, 1);

  G4Material* material_shell = CarbonFiber;

  G4LogicalVolume *logic_shell = new G4LogicalVolume( Carbon_shell_solid,
						      material_shell,
						      "single_absorber_logic",
						      0, 0, 0 );

  G4VisAttributes *visattchk2 = new G4VisAttributes();
  visattchk2->SetVisibility(true);
  visattchk2->SetForceSolid(true);
  visattchk2->SetColour(G4Colour::Black());
  logic_shell->SetVisAttributes(visattchk2);


  /* Place structural frame in logical tower volume */
  ostringstream name_shell;
  name_shell.str("");
  name_shell << _towerlogicnameprefix << "_single_absorber" << endl;

  new G4PVPlacement( 0, G4ThreeVector(0 , 0, 0),
		     logic_shell,
		     name_shell.str().c_str(),
		     single_tower_logic,
		     0, 0, overlapcheck);


  /* Place crystal in logical tower volume */
  ostringstream name_crystal;
  name_crystal.str("");
  name_crystal << _towerlogicnameprefix << "_single_crystal" << endl;

  new G4PVPlacement( 0, G4ThreeVector(0 , 0, 0),
		     logic_crystal,
		     name_crystal.str().c_str(),
		     single_tower_logic,
		     0, 0, overlapcheck);



  if ( verbosity > 0 )
    {
      cout << "PHG4CrystalCalorimeterDetector: Building logical volume for single tower done." << endl;
    }

  return single_tower_logic;
}

int
PHG4CrystalCalorimeterDetector::PlaceTower(G4LogicalVolume* eemcenvelope, G4LogicalVolume* singletower)
{

  /* Place single tower */
  ifstream istream_mapping;

  /* Open the datafile, if it won't open return an error */
  if (!istream_mapping.is_open())
    {
      istream_mapping.open( _mapping_tower_file.c_str() );
      if(!istream_mapping)
	{
	  cerr << "ERROR in PHG4CrystalCalorimeterDetector: Failed to open mapping file " << _mapping_tower_file << endl;
	  exit(1);
	}
    }

  string line_mapping;

  while ( getline( istream_mapping, line_mapping ) )
    {
      unsigned idx_j, idx_k, idx_l;
      G4double pos_x, pos_y, pos_z;
      G4double size_x, size_y, size_z;
      G4double alpha, beta, gamma;
      G4double dummy;

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
      if ( !( iss >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> alpha >> beta >> gamma >> dummy ) )
	{
	  cerr << "ERROR in PHG4CrystalCalorimeterDetector: Failed to read line in mapping file " << _mapping_tower_file << endl;
	  exit(1);
	}

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      ostringstream towername;
      towername.str("");
      towername << _towerlogicnameprefix << "_j_" << idx_j << "_k_" << idx_k;

      /* Add Geant4 units */
      pos_x = pos_x * cm;
      pos_y = pos_y * cm;
      pos_z = pos_z * cm;

      /* adjust tower position (absolute position given in mapping file) to be relative
       * to mother volume / envelope */
      pos_x -= _place_in_x;
      pos_y -= _place_in_y;
      pos_z -= _place_in_z;

      /* Place tower */
      if ( verbosity > 0 )
	{
	  cout << "PHG4CrystalCalorimeterDetector: Place tower " << towername.str() << endl;
	}

      new G4PVPlacement( 0, G4ThreeVector(pos_x, pos_y, pos_z),
			 singletower,
			 towername.str().c_str(),
			 eemcenvelope,
			 0, 0, overlapcheck);
    }

  return 0;
}
