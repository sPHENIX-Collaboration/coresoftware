/*!
 * \file PHG4MicromegasDetector.cc
 * \brief strongly inspired by code from Qinhua Huang <qinhua.huang@cea.fr>
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasDetector.h"

#include <phparameter/PHParameters.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4main/PHG4Detector.h>

#include <micromegas/CylinderGeomMicromegas.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                           // for PHNode
#include <phool/PHNodeIterator.h>                   // for PHNodeIterator
#include <phool/PHObject.h>                         // for PHObject

#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4String.hh>                       // for G4String
#include <Geant4/G4ThreeVector.hh>                  // for G4ThreeVector
#include <Geant4/G4Types.hh>                        // for G4double
#include <Geant4/G4VPhysicalVolume.hh>              // for G4VPhysicalVolume
#include <Geant4/G4VSolid.hh>                       // for G4VSolid

#include <cmath>
#include <iostream>
#include <numeric>
#include <tuple>                                    // for make_tuple, tuple
#include <utility>                                  // for pair, make_pair
#include <vector>                                   // for vector

//____________________________________________________________________________..
PHG4MicromegasDetector::PHG4MicromegasDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
{}

//_______________________________________________________________
bool PHG4MicromegasDetector::IsInDetector(G4VPhysicalVolume *volume) const
{ return m_activeVolumes.find( volume ) != m_activeVolumes.end(); }

//_______________________________________________________________
int PHG4MicromegasDetector::get_layer(G4VPhysicalVolume *volume) const
{
  const auto iter = m_activeVolumes.find( volume );
  return iter == m_activeVolumes.end() ? -1:iter->second;
}

//_______________________________________________________________
void PHG4MicromegasDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  create_materials();
  construct_micromegas(logicWorld);
  add_geometry_node();
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
  static const G4double temperature = 298.15*kelvin;
  static const G4double pressure = 1.*atmosphere;

  // FR4
  if (!G4Material::GetMaterial("mmg_FR4", false))
  {
    auto mmg_FR4 = new G4Material( "mmg_FR4", 1.860*g/cm3, 4, kStateSolid);
    mmg_FR4->AddMaterial( G4_C,  0.43550 );
    mmg_FR4->AddMaterial( G4_H,  0.03650 );
    mmg_FR4->AddMaterial( G4_O,  0.28120 );
    mmg_FR4->AddMaterial( G4_Si, 0.24680 );
  }

  // Kapton
  if (!G4Material::GetMaterial("mmg_Kapton", false))
  {
    auto mmg_Kapton = new G4Material( "mmg_Kapton", 1.420*g/cm3, 4, kStateSolid);
    mmg_Kapton->AddMaterial( G4_C, 0.6911330 );
    mmg_Kapton->AddMaterial( G4_H, 0.0263620 );
    mmg_Kapton->AddMaterial( G4_N, 0.0732700 );
    mmg_Kapton->AddMaterial( G4_O, 0.2092350);
  }

  // MMgas
  if (!G4Material::GetMaterial("mmg_Gas", false))
  {
    auto mmg_Gas = new G4Material( "mmg_Gas", 0.00170335*g/cm3, 3, kStateGas, temperature, pressure);
    mmg_Gas->AddMaterial( G4_Ar, 0.900 );
    mmg_Gas->AddMaterial( G4_C,  0.0826586 );
    mmg_Gas->AddMaterial( G4_H,  0.0173414 );
  }

  // MMMesh
  if (!G4Material::GetMaterial("mmg_Mesh", false))
  {
    auto mmg_Mesh = new G4Material( "mmg_Mesh", 2.8548*g/cm3, 5, kStateSolid);
    mmg_Mesh->AddMaterial( G4_Cr, 0.1900 );
    mmg_Mesh->AddMaterial( G4_Fe, 0.6800 );
    mmg_Mesh->AddMaterial( G4_Mn, 0.0200 );
    mmg_Mesh->AddMaterial( G4_Ni, 0.1000 );
    mmg_Mesh->AddMaterial( G4_Si, 0.0100 );
  }

  // MMStrips
  if (!G4Material::GetMaterial("mmg_Strips", false))
  { new G4Material( "mmg_Strips", 5.248414*g/cm3, G4_Cu, kStateSolid); }

  // MMResistivePaste
  if (!G4Material::GetMaterial("mmg_ResistPaste", false))
  { new G4Material( "mmg_ResistPaste", 0.77906*g/cm3, G4_C, kStateSolid); }

}

//_______________________________________________________________
void PHG4MicromegasDetector::construct_micromegas(G4LogicalVolume* logicWorld)
{
  // components enumeration
  /*
  this describes all the detector onion layers for a single side
  note that the detector is two sided
  */
  enum class Component
  {
    PCB,
    CuStrips,
    KaptonStrips,
    ResistiveStrips,
    Gas1,
    Mesh,
    Gas2,
    DriftCuElectrode,
    DriftKapton,
    DriftCarbon
  };

  // layer thickness
  // numbers from M. Vandenbroucke <maxence.vandenbroucke@cea.fr>
  const std::map<Component,float> layer_thickness =
  {
    { Component::PCB, 1.*mm },
    { Component::CuStrips, 12.*micrometer },
    { Component::KaptonStrips, 50.*micrometer },
    { Component::ResistiveStrips, 20.*micrometer },
    { Component::Gas1, 120.*micrometer },
    { Component::Mesh, 18.*0.8*micrometer }, // 0.8 correction factor is to account for the mesh denstity@18/45
    { Component::Gas2, 3.*mm },
    { Component::DriftCuElectrode, 15.*micrometer },
    { Component::DriftKapton, 50.*micrometer },
    { Component::DriftCarbon, 1.*mm }
  };

  // materials
  const std::map<Component,G4Material*> layer_material =
  {
    { Component::PCB, G4Material::GetMaterial("mmg_FR4") },
    { Component::CuStrips, G4Material::GetMaterial("mmg_Strips") },
    { Component::KaptonStrips, G4Material::GetMaterial("mmg_Kapton") },
    { Component::ResistiveStrips, G4Material::GetMaterial("mmg_ResistPaste" ) },
    { Component::Gas1, G4Material::GetMaterial( "mmg_Gas" ) },
    { Component::Mesh, G4Material::GetMaterial("mmg_Mesh") },
    { Component::Gas2, G4Material::GetMaterial( "mmg_Gas" ) },
    { Component::DriftCuElectrode, G4Material::GetMaterial("G4_Cu") },
    { Component::DriftKapton, G4Material::GetMaterial("mmg_Kapton") },
    { Component::DriftCarbon, G4Material::GetMaterial("G4_C") }
  };

  // color
  const std::map<Component, G4Colour> layer_color =
  {
    { Component::PCB, G4Colour::Green()},
    { Component::CuStrips, G4Colour::Brown()},
    { Component::KaptonStrips, G4Colour::Brown()},
    { Component::ResistiveStrips, G4Colour::Black()},
    { Component::Gas1, G4Colour::Grey()},
    { Component::Mesh, G4Colour::White()},
    { Component::Gas2, G4Colour::Grey()},
    { Component::DriftCuElectrode, G4Colour::Brown()},
    { Component::DriftKapton, G4Colour::Brown()},
    { Component::DriftCarbon, G4Colour(150/255., 75/255., 0)}
  };

  // setup layers in the correct order, going outwards from beam axis
  /* same compoment can appear multiple times. Layer names must be unique */
  using LayerDefinition = std::tuple<Component,std::string>;
  const std::vector<LayerDefinition> layer_definitions =
  {
    // inner side
    std::make_tuple( Component::DriftCarbon, "DriftCarbon_inner" ),
    std::make_tuple( Component::DriftKapton, "DriftKapton_inner" ),
    std::make_tuple( Component::DriftCuElectrode, "DriftCuElectrode_inner" ),
    std::make_tuple( Component::Gas2, "Gas2_inner" ),
    std::make_tuple( Component::Mesh, "Mesh_inner" ),
    std::make_tuple( Component::Gas1, "Gas1_inner" ),
    std::make_tuple( Component::ResistiveStrips, "ResistiveStrips_inner" ),
    std::make_tuple( Component::KaptonStrips, "KaptonStrips_inner" ),
    std::make_tuple( Component::CuStrips, "CuStrips_inner"  ),
    
    // PCB
    std::make_tuple( Component::PCB, "PCB" ),

    // outer side (= inner side, mirrored)
    std::make_tuple( Component::CuStrips, "CuStrips_outer" ),
    std::make_tuple( Component::KaptonStrips, "KaptonStrips_outer" ),
    std::make_tuple( Component::ResistiveStrips, "ResistiveStrips_outer" ),
    std::make_tuple( Component::Gas1, "Gas1_outer" ),
    std::make_tuple( Component::Mesh, "Mesh_outer" ),
    std::make_tuple( Component::Gas2, "Gas2_outer" ),
    std::make_tuple( Component::DriftCuElectrode, "DriftCuElectrode_outer" ),
    std::make_tuple( Component::DriftKapton, "DriftKapton_outer" ),
    std::make_tuple( Component::DriftCarbon, "DriftCarbon_outer" )
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
  auto cylinder_logic = new G4LogicalVolume( cylinder_solid, G4Material::GetMaterial("G4_AIR"), G4String(GetName()) );
  auto vis = new G4VisAttributes(G4Color(G4Colour::Grey()));
  vis->SetForceSolid(true);
  vis->SetVisibility(false);
  cylinder_logic->SetVisAttributes(vis);

  // add placement
  new G4PVPlacement( nullptr, G4ThreeVector(0,0,0), cylinder_logic, G4String(GetName()), logicWorld, false, 0, OverlapCheck() );

  // keep track of current layer
  int layer_index = m_FirstLayer;

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
    if( type == Component::Gas2 ) m_activeVolumes.insert( std::make_pair( component_phys, layer_index++ ) );
    else m_passiveVolumes.insert( component_phys );

    // update radius
    current_radius += thickness;
  }

  // print physical layers
  std::cout << "PHG4MicromegasDetector::ConstructMe - first layer: " << m_FirstLayer << std::endl;
  for( const auto& pair:m_activeVolumes )
  {  std::cout << "PHG4MicromegasDetector::ConstructMe - layer: " << pair.second << " volume: " << pair.first->GetName() << std::endl; }

  return;
}

//_______________________________________________________________
void PHG4MicromegasDetector::add_geometry_node()
{
  // do nothing if detector is inactive
  if( !m_Params->get_int_param("active")) return;

  // find or create geometry node
  std::string geonode_name = std::string( "CYLINDERGEOM_" ) + m_SuperDetector;
  auto geonode = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode_name);
  if (!geonode)
  {
    geonode = new PHG4CylinderGeomContainer();
    PHNodeIterator iter(topNode());
    auto runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
    auto newNode = new PHIODataNode<PHObject>(geonode, geonode_name, "PHObject");
    runNode->addNode(newNode);
  }

  // add cylinder objects
  /* one cylinder is added per physical volume. The dimention correspond to the drift volume */
  for( const auto& pair:m_activeVolumes )
  {
    // store layer and volume
    const int layer = pair.second;
    const G4VPhysicalVolume* volume_phys = pair.first;

    // get solid volume, cast to a tube
    const auto tub = dynamic_cast<const G4Tubs*>( volume_phys->GetLogicalVolume()->GetSolid() );

    // create cylinder and match geometry
    /* note: cylinder segmentation type and pitch is set in PHG4MicromegasHitReco */
    auto cylinder = new CylinderGeomMicromegas(layer);
    cylinder->set_radius( (tub->GetInnerRadius()/cm + tub->GetOuterRadius()/cm)/2 );
    cylinder->set_thickness( tub->GetOuterRadius()/cm - tub->GetInnerRadius()/cm );
    cylinder->set_zmin( -tub->GetZHalfLength()/cm );
    cylinder->set_zmax( tub->GetZHalfLength()/cm );
    geonode->AddLayerGeom(layer, cylinder);
  }

}
