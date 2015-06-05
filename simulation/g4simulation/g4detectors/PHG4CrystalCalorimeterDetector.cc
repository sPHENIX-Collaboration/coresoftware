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

//static double no_overlap = 0.00015 * cm; // added safety margin against overlaps by using same boundary between volumes


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
  _dx_front(50.19*mm),
  _dy_front(50.19*mm),
  _dx_back(59.3154545455*mm),
  _dy_back(59.3154545455*mm),
  _dz_crystal(90.000*mm),
  _materialCrystal( "G4_PbWO4" ),
  _active(1),
  _crystallogicnameprefix("eEcalCrystal"),
  _superdetector("NONE"),
  _inputFile( "/direct/phenix+u/jlab/github/sPHENIX-Fork/calibrations/CrystalCalorimeter/mapping/crystals_v002.txt" )
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
  bool isinabsorber = false;

  if (volume->GetName().find(_crystallogicnameprefix) != string::npos)
    {
      return 1;
    }
  else if ( isinabsorber )
    {
      return -1;
    }

  return 0;

}


//_______________________________________________________________________
void
PHG4CrystalCalorimeterDetector::Construct( G4LogicalVolume* logicWorld )
{
  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* ecal_envelope_cone = new G4Cons("eEcal_envelope_solid",
					    _rMin1, _rMax1,
					    _rMin2, _rMax2,
					    _dZ/2.,
					    _sPhi, _dPhi );

  G4LogicalVolume* ecal_envelope_log =  new G4LogicalVolume(ecal_envelope_cone, Air, G4String("eEcal_envelope"), 0, 0, 0);

  /* Define visualization attributes for envelope cone */
  G4VisAttributes* ecalVisAtt = new G4VisAttributes();
 // ecalVisAtt->SetVisibility(true);
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
		     ecal_envelope_log, "CrystalCalorimeter", logicWorld, 0, false, overlapcheck);

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

int
PHG4CrystalCalorimeterDetector::FillCrystalUnit(G4LogicalVolume *crystal_logic)
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
	
	//Air
	//G4Material* Air = G4Material::GetMaterial("G4_AIR");

	
	//*************************************
	//**********Build First Crystal********
	//*************************************
	
	//Crystal Dimensions determined by the _dx_front, with 0.18 subtracted out from all sides for carbon fiber
	G4double carbon_fiber_width = 0.18*mm;						//Width of the carbon fiber which surrounds the crystal
	G4double air_gap_carbon_fiber = 0.24*mm; 					//Air gap between crystal and the carbon fiber
	G4double air_gap_crystals = 0.60*mm;						//Air gap between crystal and crystal
	//Crystal Dimensions
	G4double dx_front_small = ( _dx_front - (2.0 * carbon_fiber_width ) - (2.0 * air_gap_carbon_fiber) - air_gap_crystals ) / 2.0;	
	G4double dy_front_small = ( _dy_front - (2.0 * carbon_fiber_width )  - (2.0 * air_gap_carbon_fiber) - air_gap_crystals ) / 2.0;	
	G4double dx_back_small = ( _dx_back - (2.0 * carbon_fiber_width )  - (2.0 * air_gap_carbon_fiber) - air_gap_crystals ) / 2.0;	
	G4double dy_back_small = (_dy_back - (2.0 * carbon_fiber_width )  - (2.0 * air_gap_carbon_fiber) - air_gap_crystals ) / 2.0;
	G4double dz = _dz_crystal;
	
	//cout << " | " << dx_front_small << " | " << dy_front_small << " | " << dx_back_small << " | " << dy_back_small << " | " << endl;

	//Carbon fiber dimensions
	G4double dx_front_large = dx_front_small + carbon_fiber_width;
	G4double dy_front_large = dy_front_small + carbon_fiber_width;
	G4double dx_back_large = dx_back_small + carbon_fiber_width;
	G4double dy_back_large = dy_back_small + carbon_fiber_width;
	
	//Vertices of the primary crystal
	std::vector<G4TwoVector> vertices;
	vertices.push_back(G4TwoVector( 0, 0));
	vertices.push_back(G4TwoVector( 0,  dy_front_small));
	vertices.push_back(G4TwoVector(  dx_front_small,  dy_front_small));
	vertices.push_back(G4TwoVector(  dx_front_small, 0));
 	vertices.push_back(G4TwoVector( 0, 0));
	vertices.push_back(G4TwoVector( 0,  dy_back_small));
	vertices.push_back(G4TwoVector(  dx_back_small,  dy_back_small));
	vertices.push_back(G4TwoVector(  dx_back_small, 0));
 	
	//Create the primary, irregularly shaped crystal
	G4VSolid* crystal_solid_small = new G4GenericTrap( G4String("eEcal_crystal"), 
							dz,
							vertices );
	
	G4LogicalVolume *crystal_logic_small = new G4LogicalVolume( crystal_solid_small,
								material_crystal,
								"eEcal_crystal",
								0, 0, 0);
	
	G4VisAttributes *visattchk_2 = new G4VisAttributes();
	visattchk_2->SetVisibility(true);
//	visattchk_2->SetVisibility(false);
	visattchk_2->SetForceSolid(true);
	visattchk_2->SetColour(G4Colour::Cyan());
	crystal_logic_small->SetVisAttributes(visattchk_2);
	
	//Create the carbon fiber backbone to surround the primary, irregularly shaped crystal on the sloping sides

	//Carbon fiber with same vertices as the primary crystal -> use same vertices
	G4VSolid* carbon_solid_small = new G4GenericTrap( G4String("carbon_solid_small"), 
							dz,
							vertices );	

	//Carbon fiber with larger dimensions (vertices + carbon fiber)
	std::vector<G4TwoVector> vertices_carbon_large;
	vertices_carbon_large.push_back(G4TwoVector( 0, 0));
	vertices_carbon_large.push_back(G4TwoVector( 0,  dy_front_large));
	vertices_carbon_large.push_back(G4TwoVector(  dx_front_large,  dy_front_large));
	vertices_carbon_large.push_back(G4TwoVector(  dx_front_large, 0));
 	vertices_carbon_large.push_back(G4TwoVector( 0, 0));
	vertices_carbon_large.push_back(G4TwoVector( 0,  dy_back_large));
	vertices_carbon_large.push_back(G4TwoVector(  dx_back_large,  dy_back_large));
	vertices_carbon_large.push_back(G4TwoVector(  dx_back_large, 0));
	
	G4VSolid* carbon_solid_large = new G4GenericTrap( G4String("carbon_solid_large"), 
							dz,
							vertices_carbon_large );
	
	//Use G4SubtractionSolid to subtract 1 from 2
	G4SubtractionSolid* Carbon_Shell_Solid = new G4SubtractionSolid(G4String("Carbon_Shell_Solid"), 
								  carbon_solid_large,
								  carbon_solid_small);
	
	G4LogicalVolume *Carbon_Shell = new G4LogicalVolume( Carbon_Shell_Solid,
							CarbonFiber,
							"Carbon_Shell",
							0, 0, 0);

        G4VisAttributes *visattchk_3 = new G4VisAttributes();
//        visattchk_3->SetVisibility(true);
	visattchk_3->SetVisibility(false);
        visattchk_3->SetForceSolid(true);
        visattchk_3->SetColour(G4Colour::Black());
        Carbon_Shell->SetVisAttributes(visattchk_3);

	
	//Read in mapping file for a single 4 x 4 block

	const string Crystal_Mapping_Small = "/direct/phenix+u/jlab/github/wip/Crystal_Subunit_Mapping.txt";
	const int NumberOfIndices = 6; // Number of indices in mapping file for 4x4 block

        ifstream datafile_2;

        if (!datafile_2.is_open())
        {
                datafile_2.open(Crystal_Mapping_Small.c_str());
                if(!datafile_2)
                {
                        cerr << endl << "*******************************************************************" << endl;
                        cerr << "ERROR in 4 by 4 crystal mapping"; 
			cerr << "Failed to open " << Crystal_Mapping_Small <<" --- Exiting program." << endl;
                        cerr << "*******************************************************************" << endl << endl;
                        exit(1);
                }
        }

	int NumberOfLines = 0;
        ifstream in(Crystal_Mapping_Small.c_str());
        std::string unused;
        while ( std::getline(in, unused) )
           ++NumberOfLines;

	G4int j_cry = NumberOfLines;
	G4int k_cry = NumberOfIndices;

	double FourByFour[j_cry][k_cry];

	G4int j = 0;
        G4int k = 0;

        while (j_cry > j) {
                while (k_cry > k) {
                        datafile_2 >> FourByFour[j][k];
                        k++;
                }
                j++;
                k = 0;
        }
	
	//Place this new mother volume 16 times in the mother volume (called in function)
	
	G4int j_idx, k_idx;
	G4double x_cent, y_cent, z_cent, rot_z;

	j = 0;
	while (j_cry > j) {
		j_idx = FourByFour[j][0];
		k_idx = FourByFour[j][1];
		x_cent = FourByFour[j][2];
		y_cent = FourByFour[j][3];
		z_cent = FourByFour[j][4];
		rot_z = FourByFour[j][5];

	        G4ThreeVector Crystal_Center = G4ThreeVector(x_cent*mm, y_cent*mm, z_cent*mm);

	        G4RotationMatrix *Rot = new G4RotationMatrix(); //rotation matrix for the placement of each crystal
		        Rot->rotateX(0*rad);
		        Rot->rotateY(0*rad);
		        Rot->rotateZ(rot_z*rad);

		ostringstream crystal_name;
		crystal_name.str("");
	        crystal_name << "eEcal_crystal" << "_j_"<< j_idx << "_k_" << k_idx;

		ostringstream carbon_fiber_name;
		carbon_fiber_name.str("");
		carbon_fiber_name << "Carbon_Fiber_" << j;

		new G4PVPlacement( Rot, Crystal_Center,
        	        crystal_logic_small,
        	        crystal_name.str().c_str(),
        	        crystal_logic,
        	        0, 0, overlapcheck);

                G4ThreeVector Carbon_Center = G4ThreeVector((x_cent + carbon_fiber_width / 2.0 )*mm, (y_cent + carbon_fiber_width / 2.0 )*mm, z_cent*mm); //not sure about this placement.. have to check visualization.
		new G4PVPlacement(Rot, Carbon_Center,
			Carbon_Shell,
			carbon_fiber_name.str().c_str(),
			crystal_logic,
			0, 0, overlapcheck);

		j_idx = k_idx = 0;
		x_cent = y_cent = z_cent = rot_z = 0.0;
		j++;
	}
	
	return 0;
	

}

//_______________________________________________________________________
int
PHG4CrystalCalorimeterDetector::ConstructCrystals(G4LogicalVolume* ecalenvelope)
{

	G4int NumberOfLines;		 			//Number of crystals to be created.
	const G4int NumberOfIndices = 7; 			//Different dimensions needed for crystal placement
	const string FileName = _inputFile.c_str();		//File in which crystal positions are stored
	
	G4int j_cry, k_cry;					//Indices for matrix
	G4int j_idx, k_idx;					//Indices of each crstals
	G4double x_cent, y_cent, z_cent, r_theta, r_phi;	//Coordinates of crystal in [x,y,z,theta,phi]

	G4double dx1 = 0.0;					//Half of the extent of the front face of the trapezoid in x
	G4double dx2 = 0.0; 					//Half of the extent of the back face of the trapezoid in x
	G4double dy1 = 0.0;					//Half of the extent of the front face of the trapezoid in y
	G4double dy2 = 0.0;					//Half of the extent of the back face of the trapezoid in y
	G4double dz =  0.0;					//Half of the extent of the crystal in z

	CrystalDimensions(dx1, dy1, dx2, dy2, dz); 		//File crystal dimensions with function PHG4CrystalCalorimeterDetector::CrystalDimensions

	//Create single crystal

	G4VSolid* crystal_solid = new G4Trd(G4String("eEcal_crystal"),
		dx1,						//Half length on the small face in x
		dx2,						//Half length on the large face in x
		dy1,						//Half length on the small face in y
		dy2,						//Half length on the large face in y
		dz);						//Half length in z

//	G4Material* material_crystal = G4Material::GetMaterial(_materialCrystal.c_str());
	G4Material* Air = G4Material::GetMaterial("G4_AIR");

	G4LogicalVolume *crystal_logic = new G4LogicalVolume( crystal_solid,
		Air,
		"eEcal_crystal_unit",
		0, 0, 0);

	G4VisAttributes *visattchk = new G4VisAttributes();
//	visattchk->SetVisibility(true);
	visattchk->SetVisibility(false);
	visattchk->SetForceSolid(true);
	visattchk->SetColour(G4Colour::Yellow());
	crystal_logic->SetVisAttributes(visattchk);
	
	FillCrystalUnit(crystal_logic);

	ostringstream name;

	ifstream datafile;

	if (!datafile.is_open())
	{
		datafile.open(FileName.c_str());
		if(!datafile)
		{
			cerr << endl << "*******************************************************************" << endl;
			cerr << "ERROR: Failed to open " << FileName <<" --- Exiting program." << endl;
			cerr << "*******************************************************************" << endl << endl;
			exit(1);
		}
	}

	//Determine the number of crystals to be created
	NumberOfLines = 0;
	ifstream in(_inputFile.c_str());
	std::string unused;
	while ( std::getline(in, unused) )
	   ++NumberOfLines;

	j_cry = NumberOfLines; 		// = Number of Crystals
	k_cry = NumberOfIndices; 	// = j, k, x, y, z, alpha, beta.
	
	double Crystals[j_cry][k_cry];

	G4int j = 0;
	G4int k = 0;

	//Fill matrix with the data from the mapping file
	while (j_cry > j) {
		while (k_cry > k) {
			datafile >> Crystals[j][k];
			k++;
		}
		j++;
		k = 0;
	}

	// Find the maximum k index
	G4int k_max = 0;
	j = 0;
	while (j_cry > j){
		if(Crystals[j][1] > k_max) k_max = Crystals[j][1];
		j++;
	}
	
	//Second Quadrant
	j = 0;
	while (j_cry > j) {
		j_idx = Crystals[j][0];
		k_idx = Crystals[j][1];
		x_cent = Crystals[j][2] - _place_in_x;
		y_cent = Crystals[j][3] - _place_in_y;
		z_cent = Crystals[j][4] - _place_in_z; 	//Coordinate system refers to mother volume, have to subtract out its position in the actual xyz-space
		r_theta = Crystals[j][5];		//Rotation in Horizontal
		r_phi = Crystals[j][6];			//Rotation in Vertical

		G4ThreeVector Crystal_Center = G4ThreeVector(x_cent*mm, y_cent*mm, z_cent*mm);

		G4RotationMatrix *Rot = new G4RotationMatrix(); //rotation matrix for the placement of each crystal
		Rot->rotateX(r_phi*rad);
		Rot->rotateY(r_theta*rad);
		Rot->rotateZ(0*rad);

		name.str("");
		name << _crystallogicnameprefix << "_j_" << j_idx << "_k_" << k_idx;
		
		new G4PVPlacement( Rot, Crystal_Center,
			crystal_logic,
			name.str().c_str(),
			ecalenvelope,
			0, 0, overlapcheck);

		j_idx = k_idx = 0;
		x_cent = y_cent = z_cent = r_theta = r_phi = 0.0;

		j++;
	}


	//First Quadrant
        j = 0;
        while (j_cry > j) {
                j_idx = k_max - Crystals[j][0];
                k_idx = Crystals[j][1];
                x_cent = -1.0 * ( Crystals[j][2] - _place_in_x );
                y_cent = Crystals[j][3] - _place_in_y;
                z_cent = Crystals[j][4] - _place_in_z;
                r_theta = -1.0*Crystals[j][5];
                r_phi = Crystals[j][6];

                G4ThreeVector Crystal_Center = G4ThreeVector(x_cent*mm, y_cent*mm, z_cent*mm);

                G4RotationMatrix *Rot = new G4RotationMatrix(); //rotation matrix for the placement of each crystal
	        Rot->rotateX(r_phi*rad);
                Rot->rotateY(r_theta*rad);
                Rot->rotateZ(0*rad);

                name.str("");
                name << _crystallogicnameprefix << "_j_" << j_idx << "_k_" << k_idx;

                new G4PVPlacement( Rot, Crystal_Center,
                        crystal_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, overlapcheck);

                j_idx = k_idx = 0;
                x_cent = y_cent = z_cent = r_theta = r_phi = 0.0;

                j++;
        }

	//Fourth Quadrant
        j = 0;
        while (j_cry > j) {
                j_idx = Crystals[j][0];
                k_idx = k_max - Crystals[j][1];
                x_cent = -1.0 * ( Crystals[j][2] - _place_in_x );
                y_cent = -1.0 * ( Crystals[j][3] - _place_in_y );
                z_cent = Crystals[j][4] - _place_in_z;
                r_theta = -1.0*Crystals[j][5];
                r_phi = -1.0*Crystals[j][6];

                G4ThreeVector Crystal_Center = G4ThreeVector(x_cent*mm, y_cent*mm, z_cent*mm);

		G4RotationMatrix *Rot = new G4RotationMatrix(); //rotation matrix for the placement of each crystal
                Rot->rotateX(r_phi*rad);
                Rot->rotateY(r_theta*rad);
                Rot->rotateZ(0*rad);

                name.str("");
                name << _crystallogicnameprefix << "_j_" << j_idx << "_k_" << k_idx;

                new G4PVPlacement( Rot, Crystal_Center,
                        crystal_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, overlapcheck);

                j++;
        }

	//Third Quadrant
        j = 0;
        while (j_cry > j) {
                j_idx = k_max - Crystals[j][0];
                k_idx = k_max - Crystals[j][1];
                x_cent = Crystals[j][2] - _place_in_x;
                y_cent = -1.0 * ( Crystals[j][3] - _place_in_y );
                z_cent = Crystals[j][4] - _place_in_z;
                r_theta = Crystals[j][5];
                r_phi = -1.0*Crystals[j][6];

                G4ThreeVector Crystal_Center = G4ThreeVector(x_cent*mm, y_cent*mm, z_cent*mm);

		G4RotationMatrix *Rot = new G4RotationMatrix(); //rotation matrix for the placement of each crystal
                Rot->rotateX(r_phi*rad);
                Rot->rotateY(r_theta*rad);
                Rot->rotateZ(0*rad);

                name.str("");
                name << _crystallogicnameprefix << "_j_" << j_idx << "_k_" << k_idx;

                new G4PVPlacement( Rot, Crystal_Center,
                        crystal_logic,
                        name.str().c_str(),
                        ecalenvelope,
                        0, 0, overlapcheck);

                j++;
        }
	
	return 0;
}
