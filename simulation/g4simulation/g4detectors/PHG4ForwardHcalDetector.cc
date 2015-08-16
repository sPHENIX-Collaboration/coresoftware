#include "PHG4ForwardHcalDetector.h"
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
PHG4ForwardHcalDetector::PHG4ForwardHcalDetector( PHCompositeNode *Node, const std::string &dnam ):
  PHG4Detector(Node, dnam),
  _place_in_x(0.0*mm),
  _place_in_y(0.0*mm),
  _place_in_z(4000.0*mm),
  _rot_in_x(0.0),
  _rot_in_y(0.0),
  _rot_in_z(0.0),
  _rMin1(50*mm),
  _rMax1(2620*mm),
  _rMin2(50*mm),
  _rMax2(3369*mm),
  _dZ(1000*mm),
  _sPhi(0),
  _dPhi(2*M_PI),
  _tower_dx(100*mm),
  _tower_dy(100*mm),
  _tower_dz(1000.0*mm),
  _materialScintillator( "G4_PLASTIC_SC_VINYLTOLUENE" ),
  _materialAbsorber( "G4_Fe" ),
  _active(1),
  _towerlogicnameprefix("hHcalTower"),
  _superdetector("NONE"),
  _mapping_tower_file("")
{

}


//_______________________________________________________________________
PHG4ForwardHcalDetector::~PHG4ForwardHcalDetector()
{}


//_______________________________________________________________________
int
PHG4ForwardHcalDetector::IsInForwardHcal(G4VPhysicalVolume * volume) const
{
  if (volume->GetName().find(_towerlogicnameprefix) != string::npos)
    {
      if (volume->GetName().find("scintillator") != string::npos)
	{
	  return 1;
	}
      /* only record energy in actual absorber- drop energy lost in air gaps inside hcal envelope */
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
PHG4ForwardHcalDetector::Construct( G4LogicalVolume* logicWorld )
{
  if ( verbosity > 0 )
    {
      cout << "PHG4ForwardHcalDetector: Begin Construction" << endl;
    }


  if ( _mapping_tower_file.empty() )
    {
      cout << "ERROR in PHG4ForwardHcalDetector: No tower mapping file specified. Abort detector construction." << endl;
      cout << "Please run SetTowerMappingFile( std::string filename ) first." << endl;
      exit(1);
    }

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* hcal_envelope_solid = new G4Cons("hHcal_envelope_solid",
					    _rMin1, _rMax1,
					    _rMin2, _rMax2,
					    _dZ/2.,
					    _sPhi, _dPhi );

  G4LogicalVolume* hcal_envelope_log =  new G4LogicalVolume(hcal_envelope_solid, Air, G4String("hHcal_envelope"), 0, 0, 0);

  /* Define visualization attributes for envelope cone */
  G4VisAttributes* hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(false);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Magenta());
  hcal_envelope_log->SetVisAttributes(hcalVisAtt);

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(_rot_in_x);
  hcal_rotm.rotateY(_rot_in_y);
  hcal_rotm.rotateZ(_rot_in_z);

  /* Place envelope cone in simulation */
  ostringstream name_envelope;
  name_envelope.str("");
  name_envelope << _towerlogicnameprefix << "_envelope" << endl;

  new G4PVPlacement( G4Transform3D(hcal_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z) ),
		     hcal_envelope_log, name_envelope.str().c_str(), logicWorld, 0, false, overlapcheck);

  /* Construct single calorimeter tower */
  G4LogicalVolume* singletower = ConstructTower();

  /* Place calorimeter tower within envelope */
  PlaceTower( hcal_envelope_log , singletower );

  return;
}


//_______________________________________________________________________
G4LogicalVolume*
PHG4ForwardHcalDetector::ConstructTower()
{
  if ( verbosity > 0 )
    {
      cout << "PHG4ForwardHcalDetector: Build logical volume for single tower..." << endl;
    }

  /* create logical volume for single tower */
  G4Material* material_air = G4Material::GetMaterial( "G4_AIR" );

  G4VSolid* single_tower_solid = new G4Box( G4String("single_tower_solid"),
                                       _tower_dx / 2.0,
                                       _tower_dy / 2.0,
                                       _tower_dz / 2.0 );

  G4LogicalVolume *single_tower_logic = new G4LogicalVolume( single_tower_solid,
						      material_air,
						      "single_tower_logic",
						      0, 0, 0);

  /* create geometry volumes for scintillator and absorber plates to place inside single_tower */
  G4int nlayers = 30;
  G4double thickness_layer = _tower_dz/(float)nlayers;
  G4double thickness_absorber = thickness_layer / 5.0 * 4.0; // 4/5th absorber
  G4double thickness_scintillator = thickness_layer / 5.0 * 1.0; // 1/5th scintillator

  G4VSolid* solid_absorber = new G4Box( G4String("single_plate_absorber_solid"),
				   _tower_dx / 2.0,
				   _tower_dy / 2.0,
				   thickness_absorber / 2.0 );

  G4VSolid* solid_scintillator = new G4Box( G4String("single_plate_scintillator"),
				   _tower_dx / 2.0,
				   _tower_dy / 2.0,
				   thickness_scintillator / 2.0 );

  /* create logical volumes for scintillator and absorber plates to place inside single_tower */
  G4Material* material_scintillator = G4Material::GetMaterial( _materialScintillator.c_str() );
  G4Material* material_absorber = G4Material::GetMaterial( _materialAbsorber.c_str() );

  G4LogicalVolume *logic_absorber = new G4LogicalVolume( solid_absorber,
						    material_absorber,
						    "single_plate_absorber_logic",
						    0, 0, 0);

  G4LogicalVolume *logic_scint = new G4LogicalVolume( solid_scintillator,
						    material_scintillator,
						    "hHcal_scintillator_plate_logic",
						    0, 0, 0);


  /* place physical volumes for absorber and scintillator plates */
  G4double xpos_i = 0;
  G4double ypos_i = 0;
  G4double zpos_i = ( -1 * _tower_dz / 2.0 ) + thickness_absorber / 2.0;

  ostringstream name_absorber;
  name_absorber.str("");
  name_absorber << _towerlogicnameprefix << "_single_plate_absorber" << endl;

  ostringstream name_scintillator;
  name_scintillator.str("");
  name_scintillator << _towerlogicnameprefix << "_single_plate_scintillator" << endl;

  for (int i = 1; i <= nlayers; i++)
    {
      new G4PVPlacement( 0, G4ThreeVector(xpos_i , ypos_i, zpos_i),
			 logic_absorber,
			 name_absorber.str().c_str(),
			 single_tower_logic,
			 0, 0, overlapcheck);

      zpos_i += ( thickness_absorber/2. + thickness_scintillator/2. );

      new G4PVPlacement( 0, G4ThreeVector(xpos_i , ypos_i, zpos_i),
			 logic_scint,
			 name_scintillator.str().c_str(),
			 single_tower_logic,
			 0, 0, overlapcheck);

      zpos_i += ( thickness_absorber/2. + thickness_scintillator/2. );
    }


  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(true);
  visattchk->SetColour(G4Colour::Cyan());
  single_tower_logic->SetVisAttributes(visattchk);

  if ( verbosity > 0 )
    {
      cout << "PHG4ForwardHcalDetector: Building logical volume for single tower done." << endl;
    }

  return single_tower_logic;
}

int
PHG4ForwardHcalDetector::PlaceTower(G4LogicalVolume* hcalenvelope, G4LogicalVolume* singletower)
{

  /* Place single tower */
  ifstream istream_mapping;

  /* Open the datafile, if it won't open return an error */
  if (!istream_mapping.is_open())
    {
      istream_mapping.open( _mapping_tower_file.c_str() );
      if(!istream_mapping)
	{
	  cerr << "ERROR in PHG4ForwardHcalDetector: Failed to open mapping file " << _mapping_tower_file << endl;
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
	      cout << "PHG4ForwardHcalDetector: SKIPPING line in mapping file: " << line_mapping << endl;
	    }
	  continue;
	}

      /* read string- break if error */
      if ( !( iss >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> alpha >> beta >> gamma >> dummy ) )
	{
	  cerr << "ERROR in PHG4ForwardHcalDetector: Failed to read line in mapping file " << _mapping_tower_file << endl;
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
	  cout << "PHG4ForwardHcalDetector: Place tower " << towername.str() << endl;
	}

      new G4PVPlacement( 0, G4ThreeVector(pos_x , pos_y, pos_z),
			 singletower,
			 towername.str().c_str(),
			 hcalenvelope,
			 0, 0, overlapcheck);
    }

  return 0;
}
