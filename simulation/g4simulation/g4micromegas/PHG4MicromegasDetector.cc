/*!
 * \file PHG4MicromegasDetector.cc
 * \brief strongly inspired by code from Qinhua Huang <qinhua.huang@cea.fr>
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

//____________________________________________________________________________..
//
// This is a working template for the G4 Construct() method which needs to be implemented
// We wedge a method between the G4 Construct() to enable volume hierarchies on the macro
// so here it is called ConstructMe() but there is no functional difference
// Currently this installs a simple G4Box solid, creates a logical volume from it
// and places it. Put your own detector in place (just make sure all active volumes
// get inserted into the m_PhysicalVolumes)
//
// Rather than using hardcoded values you should consider using the parameter class
// Parameter names and defaults are set in PHG4MicromegasSubsystem::SetDefaultParameters()
// Only parameters defined there can be used (also to override in the macro)
// to avoids typos.
// IMPORTANT: parameters have no inherent units, there is a convention (cm/deg)
// but in any case you need to multiply them here with the correct CLHEP/G4 unit
//
// The place where you put your own detector is marked with
// //begin implement your own here://
// //end implement your own here://
// Do not forget to include the G4 includes for your volumes
//____________________________________________________________________________..

#include "PHG4MicromegasDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4Subsystem.h>
#include <phool/phool.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <numeric>

class G4VSolid;
class PHCompositeNode;

//____________________________________________________________________________..
PHG4MicromegasDetector::PHG4MicromegasDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
{}

//_______________________________________________________________
bool PHG4MicromegasDetector::IsInDetector(G4VPhysicalVolume *volume) const
{ return m_PhysicalVolumes.find( volume ) != m_PhysicalVolumes.end(); }

//_______________________________________________________________
int PHG4MicromegasDetector::get_layer(G4VPhysicalVolume *volume) const
{
  const auto iter = m_PhysicalVolumes.find( volume );
  return iter == m_PhysicalVolumes.end() ? -1:iter->second;
}

//_______________________________________________________________
void PHG4MicromegasDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  create_materials();

  // keep
  auto gasMaterial = G4Material::GetMaterial( "myMMGas" );

  // components enumeration
  /*
  this describes all the detector onion layers for a single side
  note that the detector is two sided
  */
  enum Component
  {
    CuGround,
    PCB,
    CuStrips,
    KaptonStrips,
    ResistiveStrips,
    Gas1,
    Mesh,
    Gas2,
    DriftCuElectrode,
    DriftKapton,
    DriftCuGround,
  };

  // layer thickness
  const std::map<Component,float> layer_thickness =
  {
    { CuGround, 0.000158*cm },
    { PCB, 0.01*cm },
    { CuStrips, 0.0012*cm },
    { KaptonStrips, 0.0075*cm },
    { ResistiveStrips, 0.002*cm },
    { Gas1, 0.002*cm },
    { Mesh, 0.0018*cm },
    { Gas2, 0.3*cm },
    { DriftCuElectrode, 0.0005*cm },
    { DriftKapton, 0.025*cm },
    { DriftCuGround, 0.000041*cm },
  };

  // materials
  const std::map<Component,G4Material*> layer_material =
  {
    { CuGround, G4Material::GetMaterial("myCopper") },
    { PCB, G4Material::GetMaterial("myFR4") },
    { CuStrips, G4Material::GetMaterial("myMMStrips") },
    { KaptonStrips, G4Material::GetMaterial("myKapton") },
    { ResistiveStrips, G4Material::GetMaterial("myMMResistivePaste" ) },
    { Gas1, gasMaterial },
    { Mesh, G4Material::GetMaterial("myMMMesh") },
    { Gas2, gasMaterial },
    { DriftCuElectrode, G4Material::GetMaterial("myCopper") },
    { DriftKapton, G4Material::GetMaterial("myKapton") },
    { DriftCuGround, G4Material::GetMaterial("myCopper") }
  };

  // color
  const std::map<int, G4Colour> layer_color =
  {
    { CuGround, G4Colour::Brown()},
    { PCB, G4Colour::Green()},
    { CuStrips, G4Colour::Brown()},
    { KaptonStrips, G4Colour::Brown()},
    { ResistiveStrips, G4Colour::Black()},
    { Gas1, G4Colour::Grey()},
	   { Mesh, G4Colour::White()},
	   { Gas2, G4Colour::Grey()},
	   { DriftCuElectrode, G4Colour::Brown()},
	   { DriftKapton, G4Colour::Brown()},
    { DriftCuGround, G4Colour(51/255., 26/255., 0)}
  };

  // setup layers in the correct order, going outwards from beam axis
  using LayerDefinition = std::tuple<Component,std::string>;
  const std::vector<LayerDefinition> layer_definitions =
  {
    // inner side
    { DriftCuGround, "DriftCuGround_inner"},
    { DriftKapton, "DriftKapton_inner"},
    { DriftCuElectrode, "DriftCuElectrode_inner"},
    { Gas2, "Gas2_inner"},
    { Mesh, "Mesh_inner"},
    { Gas1, "Gas1_inner"},
    { ResistiveStrips, "ResistiveStrips_inner"},
    { KaptonStrips, "KaptonStrips_inner"},
    { CuStrips, "CuStrips_inner"},
    { PCB, "PCB_inner"},

    // separating ground
    { CuGround, "CuGround"},

    // outer side (= inner side, mirrored)
    { PCB, "PCB_outer"},
    { CuStrips, "CuStrips_outer"},
    { KaptonStrips, "KaptonStrips_outer"},
    { ResistiveStrips, "ResistiveStrips_outer"},
    { Gas1, "Gas1_outer"},
    { Mesh, "Mesh_outer"},
    { Gas2, "Gas2_outer"},
    { DriftCuElectrode, "DriftCuElectrode_outer"},
    { DriftKapton, "DriftKapton_outer"},
    { DriftCuGround, "DriftCuGround_outer"}
  };

  // start seting up volumes
  // get initial radius
  const double radius = m_Params->get_double_param("mm_radius")*cm;
  const double length =  m_Params->get_double_param("mm_length")*cm;

  // get total thickness
  const double thickness = std::accumulate(
    layer_definitions.begin(), layer_definitions.end(), 0.,
    [layer_thickness](double value, LayerDefinition layer )
    { return value + layer_thickness.at(std::get<0>(layer)); } );

  std::cout << "PHG4MicromegasDetector::ConstructMe - detector thickness is " << thickness/cm << " cm" << std::endl;

  // create mother volume
  auto cylinder_solid = new G4Tubs( G4String(GetName()), radius - 0.001*mm, radius + thickness + 0.001*mm, length / 2., 0, M_PI*2);
  auto cylinder_logic = new G4LogicalVolume( cylinder_solid, G4Material::GetMaterial("myAir"), G4String(GetName()) );
  auto vis = new G4VisAttributes(G4Color(G4Colour::Grey()));
  vis->SetForceSolid(true);
  vis->SetVisibility(false);
  cylinder_logic->SetVisAttributes(vis);

  // add placement
  new G4PVPlacement( nullptr, G4ThreeVector(0,0,0), cylinder_logic, G4String(GetName()), logicWorld, false, 0, OverlapCheck() );

  // keep track of current layer
  int layer_index = m_first_layer;

  // create detector
  /* we loop over registered layers and create volumes for each */
  auto current_radius = radius;
  for( const auto& layer:layer_definitions )
  {
    const Component& type = std::get<0>(layer);
    const std::string& name = std::get<1>(layer);

    // layer name
    G4String cname = G4String(GetName()) + "_" + name;

    // get thickness, material and name
    const auto thickness = layer_thickness.at(type);
    const auto material = layer_material.at(type);
    const auto color = layer_color.at(type);

    auto component_solid = new G4Tubs(cname+"_solid", current_radius, current_radius+thickness, length/2, 0, M_PI*2);
    auto component_logic = new G4LogicalVolume( component_solid, material, cname+"_logic");
    auto vis = new G4VisAttributes( color );
    vis->SetForceSolid(true);
    vis->SetVisibility(true);
    component_logic->SetVisAttributes(vis);

    auto component_phys = new G4PVPlacement( nullptr, G4ThreeVector(0,0,0), component_logic, cname+"_phys", cylinder_logic, false, 0, OverlapCheck() );

    // store active volume
    if( type == Gas2 ) m_PhysicalVolumes.insert( std::make_pair( component_phys, layer_index++ ) );

    // update radius
    current_radius += thickness;
  }

  // print physical layers
  std::cout << "PHG4MicromegasDetector::ConstructMe - first layer: " << m_first_layer << std::endl;
  for( const auto& pair:m_PhysicalVolumes )
  {  std::cout << "PHG4MicromegasDetector::ConstructMe - layer: " << pair.second << " volume: " << pair.first->GetName() << std::endl; }

  return;
}

//_______________________________________________________________
void PHG4MicromegasDetector::Print(const std::string &what) const
{
  std::cout << "PHG4Micromegas Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Version 0.1" << std::endl;
    std::cout << "Parameters:" << std::endl;
    m_Params->Print();
  }
  return;
}

//_______________________________________________________________
void PHG4MicromegasDetector::create_materials() const
{
  // get the list of NIST materials
  // ---------------------------------
  auto G4_N = G4Material::GetMaterial("G4_N");
  auto G4_O = G4Material::GetMaterial("G4_O");
  auto G4_C = G4Material::GetMaterial("G4_C");
  auto G4_H = G4Material::GetMaterial("G4_H");
  auto G4_Si = G4Material::GetMaterial("G4_Si");
  auto G4_Ar = G4Material::GetMaterial("G4_Ar");
  auto G4_Cr = G4Material::GetMaterial("G4_Cr");
  auto G4_Fe = G4Material::GetMaterial("G4_Fe");
  auto G4_Mn = G4Material::GetMaterial("G4_Mn");
  auto G4_Ni = G4Material::GetMaterial("G4_Ni");
  auto G4_Cu = G4Material::GetMaterial("G4_Cu");

  // combine elements
  // ----------------
  static constexpr G4double temperature = 298.15*kelvin;
  static constexpr G4double pressure = 1.*atmosphere;


  // air
  if (!G4Material::GetMaterial("myAir", false))
  {
    auto myAir = new G4Material( "myAir", 0.001205*g/cm3, 2, kStateGas, temperature, pressure);
    myAir->AddMaterial( G4_N, 0.77 );
    myAir->AddMaterial( G4_O, 0.23 );
  }

  // FR4
  if (!G4Material::GetMaterial("myFR4", false))
  {
    auto myFR4 = new G4Material( "myFR4", 1.860*g/cm3, 4, kStateSolid);
    myFR4->AddMaterial( G4_C,  0.43550 );
    myFR4->AddMaterial( G4_H,  0.03650 );
    myFR4->AddMaterial( G4_O,  0.28120 );
    myFR4->AddMaterial( G4_Si, 0.24680 );
  }

  // Kapton
  if (!G4Material::GetMaterial("myKapton", false))
  {
    auto myKapton = new G4Material( "myKapton", 1.420*g/cm3, 4, kStateSolid);
    myKapton->AddMaterial( G4_C, 0.6911330 );
    myKapton->AddMaterial( G4_H, 0.0263620 );
    myKapton->AddMaterial( G4_N, 0.0732700 );
    myKapton->AddMaterial( G4_O, 0.2092350);
  }

  // MMgas
  if (!G4Material::GetMaterial("myMMGas", false))
  {
    auto myMMGas = new G4Material( "myMMGas", 0.00170335*g/cm3, 3, kStateGas, temperature, pressure);
    myMMGas->AddMaterial( G4_Ar, 0.900 );
    myMMGas->AddMaterial( G4_C,  0.0826586 );
    myMMGas->AddMaterial( G4_H,  0.0173414 );
  }

  // MMMesh
  if (!G4Material::GetMaterial("myMMMesh", false))
  {
    auto myMMMesh = new G4Material( "myMMMesh", 2.8548*g/cm3, 5, kStateSolid);
    myMMMesh->AddMaterial( G4_Cr, 0.1900 );
    myMMMesh->AddMaterial( G4_Fe, 0.6800 );
    myMMMesh->AddMaterial( G4_Mn, 0.0200 );
    myMMMesh->AddMaterial( G4_Ni, 0.1000 );
    myMMMesh->AddMaterial( G4_Si, 0.0100 );
  }

  // MMStrips
  if (!G4Material::GetMaterial("myMMStrips", false))
  { new G4Material( "myMMStrips", 5.248414*g/cm3, G4_Cu, kStateSolid); }

  // MMResistivePaste
  if (!G4Material::GetMaterial("myMMResistivePaste", false))
  { new G4Material( "myMMResistivePaste", 0.77906*g/cm3, G4_C, kStateSolid); }

  // Copper
  if (!G4Material::GetMaterial("myCopper", false))
  { new G4Material("myCopper", 8.9600*g/cm3, G4_Cu, kStateSolid); }

}
