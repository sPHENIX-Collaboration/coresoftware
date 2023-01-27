/*!
 * \file PHG4MicromegasDetector.cc
 * \brief strongly inspired by code from Qinhua Huang <qinhua.huang@cea.fr>
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasDetector.h"

#include "PHG4MicromegasDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4Subsystem.h>

#include <micromegas/CylinderGeomMicromegas.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                           // for PHNode
#include <phool/PHNodeIterator.h>                   // for PHNodeIterator
#include <phool/PHObject.h>                         // for PHObject
#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4String.hh>                       // for G4String
#include <Geant4/G4ThreeVector.hh>                  // for G4ThreeVector
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>                        // for G4double
#include <Geant4/G4VPhysicalVolume.hh>              // for G4VPhysicalVolume
#include <Geant4/G4VSolid.hh>                       // for G4VSolid

#include <cmath>
#include <iostream>
#include <numeric>
#include <tuple>                                    // for make_tuple, tuple
#include <utility>                                  // for pair, make_pair
#include <vector>                                   // for vector

namespace
{
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
}


//____________________________________________________________________________..
PHG4MicromegasDetector::PHG4MicromegasDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4MicromegasDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
{ setup_tiles(); }

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
int PHG4MicromegasDetector::get_tileid(G4VPhysicalVolume *volume) const
{
  const auto iter = m_tiles_map.find( volume );
  return iter == m_tiles_map.end() ? -1:iter->second;
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
void PHG4MicromegasDetector::setup_tiles()
{
  
  // TODO: replace with more realistic description from latest engineering drawings
  m_tiles.clear();

  // tile dimensions
  /* they correspond to the total micromegas PCB size. They must match the definitions in construct_micromegas_tile */
  static constexpr double tile_length = 54.2; // cm
  static constexpr double tile_width = 31.6;  // cm

  {
    // bottom most sector 3pi/2 has 4 modules
    static constexpr double phi = 3.*M_PI/2;
    
    // tiles z position from CAD model (THREE PANELS UPDATE 1-18-21)
    for( const double& tile_z:{ -82.8, -27.6, 27.6, 82.2 } )
    { m_tiles.emplace_back(phi, tile_z, tile_width/CylinderGeomMicromegas::reference_radius, tile_length); }
  }
  
  {
    // neighbor sectors have two modules, separated by 10cm
    for( const double& phi: { 4.*M_PI/3, 5.*M_PI/3 } )
    {
      // tiles z position from CAD model (THREE PANELS UPDATE 1-18-21)
      for( const double& tile_z:{ -37.1, 37.1 } )
      { m_tiles.emplace_back(phi, tile_z, tile_width/CylinderGeomMicromegas::reference_radius, tile_length); }
    }  
  }
}

//_______________________________________________________________
void PHG4MicromegasDetector::create_materials() const
{
  // get the list of NIST materials
  // ---------------------------------
  auto G4_N = GetDetectorMaterial("G4_N");
  auto G4_O = GetDetectorMaterial("G4_O");
  auto G4_C = GetDetectorMaterial("G4_C");
  auto G4_H = GetDetectorMaterial("G4_H");
  auto G4_Si = GetDetectorMaterial("G4_Si");
  auto G4_Ar = GetDetectorMaterial("G4_Ar");
  auto G4_Cr = GetDetectorMaterial("G4_Cr");
  auto G4_Fe = GetDetectorMaterial("G4_Fe");
  auto G4_Mn = GetDetectorMaterial("G4_Mn");
  auto G4_Ni = GetDetectorMaterial("G4_Ni");
  auto G4_Cu = GetDetectorMaterial("G4_Cu");

  // combine elements
  // ----------------
  static const G4double temperature = 298.15*kelvin;
  static const G4double pressure = 1.*atmosphere;

  // FR4
  if (!GetDetectorMaterial("mmg_FR4", false))
  {
    auto mmg_FR4 = new G4Material( "mmg_FR4", 1.860*g/cm3, 4, kStateSolid);
    mmg_FR4->AddMaterial( G4_C,  0.43550 );
    mmg_FR4->AddMaterial( G4_H,  0.03650 );
    mmg_FR4->AddMaterial( G4_O,  0.28120 );
    mmg_FR4->AddMaterial( G4_Si, 0.24680 );
  }

  // Kapton
  if (!GetDetectorMaterial("mmg_Kapton", false))
  {
    auto mmg_Kapton = new G4Material( "mmg_Kapton", 1.420*g/cm3, 4, kStateSolid);
    mmg_Kapton->AddMaterial( G4_C, 0.6911330 );
    mmg_Kapton->AddMaterial( G4_H, 0.0263620 );
    mmg_Kapton->AddMaterial( G4_N, 0.0732700 );
    mmg_Kapton->AddMaterial( G4_O, 0.2092350);
  }

  // MMgas
  if (!GetDetectorMaterial("mmg_Gas", false))
  {
    auto mmg_Gas = new G4Material( "mmg_Gas", 0.00170335*g/cm3, 3, kStateGas, temperature, pressure);
    mmg_Gas->AddMaterial( G4_Ar, 0.900 );
    mmg_Gas->AddMaterial( G4_C,  0.0826586 );
    mmg_Gas->AddMaterial( G4_H,  0.0173414 );
  }

  // MMMesh
  if (!GetDetectorMaterial("mmg_Mesh", false))
  {
    auto mmg_Mesh = new G4Material( "mmg_Mesh", 2.8548*g/cm3, 5, kStateSolid);
    mmg_Mesh->AddMaterial( G4_Cr, 0.1900 );
    mmg_Mesh->AddMaterial( G4_Fe, 0.6800 );
    mmg_Mesh->AddMaterial( G4_Mn, 0.0200 );
    mmg_Mesh->AddMaterial( G4_Ni, 0.1000 );
    mmg_Mesh->AddMaterial( G4_Si, 0.0100 );
  }

  // MMStrips
  if (!GetDetectorMaterial("mmg_Strips", false))
  { new G4Material( "mmg_Strips", 5.248414*g/cm3, G4_Cu, kStateSolid); }

  // MMResistivePaste
  if (!GetDetectorMaterial("mmg_ResistPaste", false))
  { new G4Material( "mmg_ResistPaste", 0.77906*g/cm3, G4_C, kStateSolid); }

}

//_______________________________________________________________
void PHG4MicromegasDetector::construct_micromegas(G4LogicalVolume* logicWorld)
{

  // start seting up volumes
  // Micromegas detector radius
  /* it corresponds to the radial position of the innermost surface of a Micromegas module , as measured in CAD model (THREE PANELS UPDATE 1-18-21) */
  static constexpr double inner_radius = 84.203*cm;

  /* 
   * this is the radius at the center of a module.
   * it is updated when constructing the first tile, assuming that all tiles are identical
   */
  double radius = 0;
  
  // create detector
  // loop over tiles
  for( size_t tileid = 0; tileid < m_tiles.size(); ++tileid )
  {

    // get relevant tile
    const auto& tile = m_tiles[tileid];

    // create tile master volume
    auto tile_logic = construct_micromegas_tile( tileid );
    
    // get tile thickness
    if( tileid == 0 )
    {
      const double tile_thickness = static_cast<G4Box*>(tile_logic->GetSolid())->GetXHalfLength()*2;
      radius = inner_radius + tile_thickness/2;
    }

    // place tile in master volume
    const double centerZ = tile.m_centerZ*cm;
    const double centerPhi = tile.m_centerPhi;

    G4RotationMatrix rotation;
    rotation.rotateZ( centerPhi*radian );
    
    const G4ThreeVector center(
      radius*std::cos(centerPhi),
      radius*std::sin(centerPhi),
      centerZ );

    G4Transform3D transform( rotation, center );
    
    const auto tilename = GetName() + "_tile_" + std::to_string(tileid);
    new G4PVPlacement( transform, tile_logic, tilename+"_phys", logicWorld, false, 0, OverlapCheck() );
  }
  
  // adjust active volume radius to account for world placement
  for( auto&& [layer, layer_radius]:m_layer_radius ) { layer_radius += radius/cm; }
  
  // print physical layers
  if( Verbosity() )
  {
    std::cout << "PHG4MicromegasDetector::ConstructMe - first layer: " << m_FirstLayer << std::endl;
    for( const auto& pair:m_activeVolumes )
    {  std::cout << "PHG4MicromegasDetector::ConstructMe - layer: " << pair.second << " volume: " << pair.first->GetName() << std::endl; }
  }

  return;
}

//_______________________________________________________________
G4LogicalVolume* PHG4MicromegasDetector::construct_micromegas_tile( int tileid )
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
    DriftCarbon,
    FeeSupport
  };

  // layer definition
  struct LayerStruct
  {
    // constructor
    LayerStruct( float thickness, G4Material* material, const G4Colour &color, double dy, double dz, double y_offset, double z_offset ):
      m_thickness( thickness ),
      m_material( material ),
      m_color( color ),
      m_dy( dy ),
      m_dz( dz ),
      m_y_offset( y_offset ),
      m_z_offset( z_offset )
    {}

    // thickness
    float m_thickness = 0;

    // material
    G4Material* m_material = nullptr;

    // color
    G4Colour m_color = 0;
    
    // dimension along y
    double m_dy = 0;
    
    // dimension along z
    double m_dz = 0;
    
    // center offset along y
    double m_y_offset = 0;
    
    // center offset along z
    double m_z_offset = 0;
    
  };

  // define all layers
  const std::map<Component, LayerStruct> layer_map =
  {
    /* adjusted PCB thickness so that total thickness up to Gas2 is 1.6mm, consistently with CAD model */
    { Component::PCB, LayerStruct(1.384*mm, GetDetectorMaterial("mmg_FR4"), G4Colour::Green(), 316*mm, 542*mm, 0, 0 )},
    { Component::CuStrips, LayerStruct(12.*micrometer, GetDetectorMaterial("mmg_Strips"), G4Colour::Brown(), 256*mm, 512*mm, -15*mm, 0)},
    { Component::KaptonStrips, LayerStruct(50.*micrometer, GetDetectorMaterial("mmg_Kapton"), G4Colour::Brown(), 256*mm, 512*mm, -15*mm, 0)},
    { Component::ResistiveStrips, LayerStruct(20.*micrometer, GetDetectorMaterial("mmg_ResistPaste" ), G4Colour::Black(), 256*mm, 512*mm, -15*mm, 0)},
    { Component::Gas1, LayerStruct(120.*micrometer, GetDetectorMaterial( "mmg_Gas" ), G4Colour::Grey(), 256*mm, 512*mm, -15*mm, 0)},
    /* 0.8 correction factor to thickness is to account for the mesh denstity@18/45 */
    { Component::Mesh, LayerStruct(18.*0.8*micrometer,  GetDetectorMaterial("mmg_Mesh"), G4Colour::White(), 256*mm, 512*mm, -15*mm, 0)}, 
    { Component::Gas2, LayerStruct(3.*mm, GetDetectorMaterial( "mmg_Gas" ), G4Colour::Grey(), 256*mm, 512*mm, -15*mm, 0)},
    { Component::DriftCuElectrode, LayerStruct(15.*micrometer, GetDetectorMaterial("G4_Cu"), G4Colour::Brown(), 256*mm, 512*mm, -15*mm, 0)},
    { Component::DriftKapton, LayerStruct(50.*micrometer, GetDetectorMaterial("mmg_Kapton"), G4Colour::Brown(), 256*mm, 512*mm, -15*mm, 0)},
    { Component::DriftCarbon, LayerStruct(1.*mm, GetDetectorMaterial("G4_C"), G4Colour(150/255., 75/255., 0), 256*mm, 512*mm, -15*mm, 0)},
    { Component::FeeSupport, LayerStruct(2.38*mm, GetDetectorMaterial("G4_Al"), G4Colour::Grey(), 182*mm, 542*mm, -67*mm, 0)}
  };

  // setup layers in the correct order, going outwards from beam axis
  /* same compoment can appear multiple times. Layer names must be unique */
  using LayerDefinition = std::tuple<Component,std::string>;
  const std::vector<LayerDefinition> layer_stack =
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
    std::make_tuple( Component::PCB, "PCB_inner" ),

    // outer side (= inner side, mirrored)
    std::make_tuple( Component::PCB, "PCB_outer" ),
    std::make_tuple( Component::CuStrips, "CuStrips_outer" ),
    std::make_tuple( Component::KaptonStrips, "KaptonStrips_outer" ),
    std::make_tuple( Component::ResistiveStrips, "ResistiveStrips_outer" ),
    std::make_tuple( Component::Gas1, "Gas1_outer" ),
    std::make_tuple( Component::Mesh, "Mesh_outer" ),
    std::make_tuple( Component::Gas2, "Gas2_outer" ),
    std::make_tuple( Component::DriftCuElectrode, "DriftCuElectrode_outer" ),
    std::make_tuple( Component::DriftKapton, "DriftKapton_outer" ),
    std::make_tuple( Component::DriftCarbon, "DriftCarbon_outer" ),
    
    // FEE support plate
    std::make_tuple( Component::FeeSupport, "FEE_support" )
  };
    
  // create two FEE boards up front, to get their total thickness and add to master volume
  std::array<G4LogicalVolume*, 2> fee_board_logic = 
  {
    construct_fee_board(0),
    construct_fee_board(1)
  };
  
  // fee thickness
  const double fee_thickness = static_cast<G4Box*>(fee_board_logic[0]->GetSolid())->GetXHalfLength()*2;
  
  // calculate total tile thickness
  const double tile_thickness = std::accumulate(
    layer_stack.begin(), layer_stack.end(), 0.,
    [&layer_map](double value, const LayerDefinition& layer )
    { return value + layer_map.at(std::get<0>(layer)).m_thickness; } ) + fee_thickness;

  // tile dimensions match that of the PCB layer
  const double tile_dy = layer_map.at(Component::PCB).m_dy;
  const double tile_dz = layer_map.at(Component::PCB).m_dz;
  
  // get world material to define parent volume
  auto rc = recoConsts::instance();
  auto world_material = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));

  // define tile name
  const auto tilename = GetName() + "_tile_" + std::to_string(tileid);
  
  auto tile_solid = new G4Box( tilename+"_solid", tile_thickness/2, tile_dy/2, tile_dz/2 );
  auto tile_logic = new G4LogicalVolume( tile_solid, world_material, "invisible_" + tilename + "_logic");
  GetDisplayAction()->AddVolume(tile_logic,G4Colour::Grey());

  /* we loop over registered layers and create volumes for each as daughter of the tile volume */
  auto current_radius_local = -tile_thickness/2;
  for( const auto& [type, name]:layer_stack )
  {

    // layer name
    /* 
     * for the Gas2 layers, which are the active components, we use a different volume name,
     * that match the old geometry implementation. This maximizes compatibility with previous versions
     */
    const G4String cname = (type == Component::Gas2) ? 
      "micromegas_measurement_" + name:
      G4String(GetName()) + "_" + name;
    
    // get thickness, material and name
    const auto& component( layer_map.at(type) );
    const auto& thickness = component.m_thickness;
    const auto& material = component.m_material;
    const auto& color = component.m_color;
    
    const auto& dy = component.m_dy;
    const auto& dz = component.m_dz;

    const auto& y_offset = component.m_y_offset;
    const auto& z_offset = component.m_z_offset;
    
    auto component_solid = new G4Box(cname+"_solid", thickness/2, dy/2, dz/2 );
    auto component_logic = new G4LogicalVolume( component_solid, material, cname+"_logic");
    GetDisplayAction()->AddVolume(component_logic , color);
    
    const G4ThreeVector center( (current_radius_local + thickness/2), y_offset, z_offset );
    auto component_phys = new G4PVPlacement( nullptr, center, component_logic, cname+"_phys", tile_logic, false, 0, OverlapCheck() );
    
    if( type == Component::Gas2 )
    {
      
      // store active volume
      // define layer from name
      const int layer_index = (name == "Gas2_inner") ? m_FirstLayer : m_FirstLayer+1;
      m_activeVolumes.insert( std::make_pair( component_phys, layer_index ) );
      m_tiles_map.insert( std::make_pair( component_phys, tileid ) );
      
      // store radius associated to this layer
      m_layer_radius.insert( std::make_pair( layer_index, (current_radius_local + thickness/2)/cm ) );
      m_layer_thickness.insert( std::make_pair( layer_index, thickness/cm) );

    } else m_passiveVolumes.insert( component_phys );

    // update radius
    current_radius_local += thickness;
  }  
  
  // add FEE boards
  /* offsets measured from CAD drawings */
  static constexpr double fee_y_offset = ( 316./2 - 44.7 - 141.5/2 )*mm;
  static constexpr double fee_x_offset = ( 542./2 - 113.1 - 140./2 )*mm;
  new G4PVPlacement( nullptr, {current_radius_local+fee_thickness/2, -fee_y_offset, -fee_x_offset}, fee_board_logic[0], GetName() + "_fee_0_phys", tile_logic, false, 0, OverlapCheck() );  
  new G4PVPlacement( nullptr, {current_radius_local+fee_thickness/2, -fee_y_offset, fee_x_offset}, fee_board_logic[1], GetName() + "_fee_1_phys", tile_logic, false, 0, OverlapCheck() );  
  
  // return master logical volume
  return tile_logic; 
}

//_______________________________________________________________
G4LogicalVolume* PHG4MicromegasDetector::construct_fee_board( int id )
{
 
  // layer definition
  struct LayerStruct
  {
    // constructor
    LayerStruct( const std::string& name, float thickness, G4Material* material, const G4Colour &color ):
      m_name( name ),
      m_thickness( thickness ),
      m_material( material ),
      m_color( color )
    {}

    // name
    std::string m_name;
    
    // thickness
    float m_thickness = 0;

    // material
    G4Material* m_material = nullptr;

    // color
    G4Colour m_color = 0;   
  };

  
  static constexpr double inch_to_cm = 2.54;

  /* 
   * FEE board consists of FR4 PCB, a coper layer, and an aluminium layer, for cooling. 
   * FR4 and Cu layer thickness taken from TPC description (PHG4TpcEndCapSubsystem::SetDefaultParameters)
   * Al layer correspond to cooling plate. Thickness from mechanical drawings
   */
  const std::vector<LayerStruct> layer_stack = {
    LayerStruct( "fee_pcb", 0.07*inch_to_cm*cm, GetDetectorMaterial("mmg_FR4"), G4Colour::Green() ),
    LayerStruct( "fee_cu", 35e-4*10*0.8*cm, GetDetectorMaterial("G4_Cu"), G4Colour::Brown() ),
    LayerStruct( "fee_al", 0.25*inch_to_cm*cm, GetDetectorMaterial("G4_Al"), G4Colour::Grey() ) };
    
  // calculate total tile thickness
  const double fee_thickness = std::accumulate(
    layer_stack.begin(), layer_stack.end(), 0.,
    [](double value, const LayerStruct& layer )
    { return value + layer.m_thickness; } );

  // fee dimensions match CAD drawings
  const double fee_dy = 141.51*mm;
  const double fee_dz = 140*mm;
  
  // get world material to define parent volume
  auto rc = recoConsts::instance();
  auto world_material = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));

  // define tile name
  const auto feename = GetName() + "_fee_" + std::to_string(id);
  
  auto fee_solid = new G4Box( feename+"_solid", fee_thickness/2, fee_dy/2, fee_dz/2 );
  auto fee_logic = new G4LogicalVolume( fee_solid, world_material, "invisible_" + feename + "_logic");
  GetDisplayAction()->AddVolume(fee_logic,G4Colour::Grey() );

  /* we loop over registered layers and create volumes for each as daughter of the fee volume */
  auto current_radius_local = -fee_thickness/2;
  for( const auto& layer:layer_stack )
  {

    // layer name
    /* 
     * for the Gas2 layers, which are the active components, we use a different volume name,
     * that match the old geometry implementation. This maximizes compatibility with previous versions
     */
    const G4String cname = G4String(GetName()) + "_" + layer.m_name;
    
    // get thickness, material and name
    const auto& thickness = layer.m_thickness;
    const auto& material = layer.m_material;
    const auto& color = layer.m_color;
        
    auto component_solid = new G4Box(cname+"_solid", thickness/2, fee_dy/2, fee_dz/2 );
    auto component_logic = new G4LogicalVolume( component_solid, material, cname+"_logic");
    GetDisplayAction()->AddVolume(component_logic , color);
    
    const G4ThreeVector center( (current_radius_local + thickness/2), 0, 0);
    auto component_phys = new G4PVPlacement( nullptr, center, component_logic, cname+"_phys", fee_logic, false, 0, OverlapCheck() );
    
    // store as passive
    m_passiveVolumes.insert( component_phys );

    // update radius
    current_radius_local += thickness;
  }   
  
  // return master logical volume
  return fee_logic;
  
}

//_______________________________________________________________
void PHG4MicromegasDetector::add_geometry_node()
{
  // do nothing if detector is inactive
  if( !m_Params->get_int_param("active")) return;

  // find or create geometry node
  const auto geonode_name = std::string( "CYLINDERGEOM_" ) + m_SuperDetector + "_FULL";
  auto geonode = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode_name);
  if (!geonode)
  {
    geonode = new PHG4CylinderGeomContainer();
    PHNodeIterator iter(topNode());
    auto runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
    auto newNode = new PHIODataNode<PHObject>(geonode, geonode_name, "PHObject");
    runNode->addNode(newNode);
  }

  // cylinder maximal length
  static constexpr double length = 220;

  // add cylinder objects
  /* one cylinder is added per physical volume. The dimention correspond to the drift volume */
  bool is_first = true;
  for( const auto& [layer_index, radius]:m_layer_radius )
  {
    // create cylinder and match geometry
    /* note: cylinder segmentation type and pitch is set in PHG4MicromegasHitReco */
    auto cylinder = new CylinderGeomMicromegas(layer_index);
    cylinder->set_radius( radius );
    cylinder->set_thickness( m_layer_thickness.at(layer_index) );
    cylinder->set_zmin( -length/2 );
    cylinder->set_zmax( length/2 );

    // tiles
    cylinder->set_tiles( m_tiles );

    /*
     * asign segmentation type and pitch
     * assume first layer in phi, other(s) are z
     */
    cylinder->set_segmentation_type( is_first ?
      MicromegasDefs::SegmentationType::SEGMENTATION_PHI :
      MicromegasDefs::SegmentationType::SEGMENTATION_Z );

    /*
     * assign drift direction
     * assume first layer is outward, with readout plane at the top, and others are inward, with readout plane at the bottom
     * this is used to properly implement transverse diffusion in ::process_event
     */
    cylinder->set_drift_direction( is_first ?
      MicromegasDefs::DriftDirection::OUTWARD :
      MicromegasDefs::DriftDirection::INWARD );

    // pitch
    // first layer (phi) has 1mm pitch. second layer (z) has 2mm pitch
    cylinder->set_pitch( is_first ? 0.1 : 0.2 );

    // if( Verbosity() )
    { cylinder->identify( std::cout ); }

    geonode->AddLayerGeom(layer_index, cylinder);

    is_first = false;

  }

}
