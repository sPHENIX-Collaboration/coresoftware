/*!
 * \file PHG4MicromegasDetector.cc
 * \brief strongly inspired by code from Qinhua Huang <qinhua.huang@cea.fr>
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasDetector.h"

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
#include <Geant4/G4VisAttributes.hh>
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

  // hard coded TPC definitions
  static constexpr double length = 210;
  static constexpr double tile_length = 50;
  static constexpr double tile_width = 25;

  {
    // bottom most sector 3pi/2 has 4 modules
    static constexpr double phi = 3.*M_PI/2;
    static constexpr int ntiles = 4;
    for (int i = 0; i < ntiles; ++i)
    { m_tiles.emplace_back(phi, length * ((0.5 + i) / ntiles - 0.5), tile_width/CylinderGeomMicromegas::reference_radius, tile_length); }
  }
  
  {
    // neighbor sectors have two modules, separated by 10cm
    static constexpr double zoffset = 10;
    for( const double& phi: { 4.*M_PI/3, 5.*M_PI/3 } )
    {
      m_tiles.emplace_back( phi, length*(1.5/4-0.5) - zoffset, tile_width/CylinderGeomMicromegas::reference_radius, tile_length );
      m_tiles.emplace_back( phi, length*(2.5/4-0.5) + zoffset, tile_width/CylinderGeomMicromegas::reference_radius, tile_length );
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
  // get initial radius
  const double radius = m_Params->get_double_param("mm_radius")*cm;

  // create detector
  // loop over tiles
  for( size_t tileid = 0; tileid < m_tiles.size(); ++tileid )
  {

    // get relevant tile
    const auto& tile = m_tiles[tileid];

    // get tile's dimension and orientation
    const double dy = CylinderGeomMicromegas::reference_radius*tile.m_sizePhi*cm;
    const double dz = tile.m_sizeZ*cm;

    // create tile master volume
    auto tile_logic = construct_micromegas_tile( tileid, dy, dz );
    
    // palce tile in master volume
    const double centerZ = tile.m_centerZ*cm;
    const double centerPhi = tile.m_centerPhi;
 
    // place
    /* not completely sure why one must rotate with oposite angle as that use for the translation */
    auto rotation = new G4RotationMatrix;
    rotation->rotateZ( -centerPhi*radian );

    const G4ThreeVector center(
      radius*std::cos(centerPhi),
      radius*std::sin(centerPhi),
      centerZ );
    
    const auto tilename = GetName() + "_tile_" + std::to_string(tileid);
    new G4PVPlacement( rotation, center, tile_logic, tilename+"_phys", logicWorld, false, 0, OverlapCheck() );
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
G4LogicalVolume* PHG4MicromegasDetector::construct_micromegas_tile( int tileid, double dy, double dz )
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

  // layer definition
  struct LayerStruct
  {
    // constructor
    LayerStruct( float thickness, G4Material* material, G4Colour color ):
      m_thickness( thickness ),
      m_material( material ),
      m_color( color )
    {}

    // thickness
    float m_thickness = 0;

    // material
    G4Material* m_material = nullptr;

    // color
    G4Colour m_color = 0;
  };

  // define all layers
  const std::map<Component, LayerStruct> layer_map =
  {
    { Component::PCB, {1.*mm, GetDetectorMaterial("mmg_FR4"), G4Colour::Green() }},
    { Component::CuStrips, { 12.*micrometer, GetDetectorMaterial("mmg_Strips"), G4Colour::Brown() }},
    { Component::KaptonStrips, { 50.*micrometer, GetDetectorMaterial("mmg_Kapton"), G4Colour::Brown() }},
    { Component::ResistiveStrips, { 20.*micrometer, GetDetectorMaterial("mmg_ResistPaste" ), G4Colour::Black() }},
    { Component::Gas1, { 120.*micrometer, GetDetectorMaterial( "mmg_Gas" ), G4Colour::Grey() }},
    { Component::Mesh, { 18.*0.8*micrometer,  GetDetectorMaterial("mmg_Mesh"), G4Colour::White()} }, // 0.8 correction factor to thickness is to account for the mesh denstity@18/45
    { Component::Gas2, { 3.*mm, GetDetectorMaterial( "mmg_Gas" ), G4Colour::Grey()}},
    { Component::DriftCuElectrode, { 15.*micrometer, GetDetectorMaterial("G4_Cu"), G4Colour::Brown() }},
    { Component::DriftKapton, { 50.*micrometer, GetDetectorMaterial("mmg_Kapton"), G4Colour::Brown() }},
    { Component::DriftCarbon, { 1.*mm, GetDetectorMaterial("G4_C"), G4Colour(150/255., 75/255., 0) }}
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
    
  // calculate total tile thickness
  const double tile_thickness = std::accumulate(
    layer_stack.begin(), layer_stack.end(), 0.,
    [&layer_map](double value, LayerDefinition layer )
    { return value + layer_map.at(std::get<0>(layer)).m_thickness; } );

  // get world material to define parent volume
  auto rc = recoConsts::instance();
  auto world_material = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));

  // define tile name
  const auto tilename = GetName() + "_tile_" + std::to_string(tileid);
  
  auto tile_solid = new G4Box( tilename+"_solid", tile_thickness/2, dy/2, dz/2 );
  auto tile_logic = new G4LogicalVolume( tile_solid, world_material, tilename+"_logic");
  auto vis = new G4VisAttributes(G4Colour::Grey());
  vis->SetForceSolid(true);
  vis->SetVisibility(false);

  tile_logic->SetVisAttributes(vis);

  /* we loop over registered layers and create volumes for each as daughter of the tile volume */
  auto current_radius_local = -tile_thickness/2;
  for( const auto& [type, name]:layer_stack )
  {

    // layer name
    const G4String cname = G4String(GetName()) + "_" + name;
    
    // get thickness, material and name
    const auto& thickness = layer_map.at(type).m_thickness;
    const auto& material = layer_map.at(type).m_material;
    const auto& color = layer_map.at(type).m_color;

    auto component_solid = new G4Box(cname+"_solid", thickness/2, dy/2, dz/2 );
    auto component_logic = new G4LogicalVolume( component_solid, material, cname+"_logic");
    auto vis = new G4VisAttributes( color );
    vis->SetForceSolid(true);
    vis->SetVisibility(true);
    component_logic->SetVisAttributes(vis);
    
    const G4ThreeVector center( (current_radius_local + thickness/2), 0, 0 );
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
  
  // return master logical volume
  return tile_logic;
  
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
  const double length =  m_Params->get_double_param("mm_length");

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
    /* they correspond to 256 channels along the phi direction, and 256 along the z direction, assuming 25x50 tiles */
    /* todo: calculate from tile dimensions */
    cylinder->set_pitch( is_first ? 25./256 : 50./256 );

    // if( Verbosity() )
    { cylinder->identify( std::cout ); }

    geonode->AddLayerGeom(layer_index, cylinder);

    is_first = false;

  }

}
