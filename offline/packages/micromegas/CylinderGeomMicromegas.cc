/*!
 * \file CylinderGeomMicromegas.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "CylinderGeomMicromegas.h"

#include <g4main/PHG4Hit.h>

#include <TVector3.h>

#include <cassert>

namespace
{
  // convenient square function
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  //! radius
  template<class T>
    inline constexpr T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  // bind angle to [-M_PI,+M_PI[. This is useful to avoid edge effects when making the difference between two angles
  template<class T>
    inline T bind_angle( const T& angle )
  {
    if( angle >= M_PI ) return angle - 2*M_PI;
    else if( angle < -M_PI ) return angle + 2*M_PI;
    else return angle;
  }

}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_local_from_world_coords( uint tileid, const TVector3& world_coordinates ) const
{
  assert( tileid < m_tiles.size() );

  // store world coordinates in array
  std::array<double,3> master;
  world_coordinates.GetXYZ( &master[0] );

  // convert to local coordinate
  std::array<double,3> local;
  transformation_matrix(tileid).MasterToLocal( &master[0], &local[0] );

  return TVector3( &local[0] );

}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_world_from_local_coords( uint tileid, const TVector3& local_coordinates ) const
{
  assert( tileid < m_tiles.size() );

  // store world coordinates in array
  std::array<double,3> local;
  local_coordinates.GetXYZ( &local[0] );

  // convert to local coordinate
  std::array<double,3> master;
  transformation_matrix(tileid).LocalToMaster( &local[0], &master[0] );

  return TVector3( &master[0] );

}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_world_from_local_vect( uint tileid, const TVector3& local_vect ) const
{
  assert( tileid < m_tiles.size() );

  // store world vect in array
  std::array<double,3> local;
  local_vect.GetXYZ( &local[0] );

  // convert to local coordinate
  std::array<double,3> master;
  transformation_matrix(tileid).LocalToMasterVect( &local[0], &master[0] );

  return TVector3( &master[0] );

}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_local_from_world_vect( uint tileid, const TVector3& world_vect ) const
{
  assert( tileid < m_tiles.size() );

  // store world vect in array
  std::array<double,3> master;
  world_vect.GetXYZ( &master[0] );

  // convert to local coordinate
  std::array<double,3> local;
  transformation_matrix(tileid).MasterToLocalVect( &master[0], &local[0] );

  return TVector3( &local[0] );

}

//________________________________________________________________________________
int CylinderGeomMicromegas::find_tile_cylindrical( const TVector3& world_coordinates ) const
{
  // check radius
  if( !check_radius(world_coordinates) ) return -1;

  // convert to polar coordinates
  const auto phi = std::atan2( world_coordinates.y(), world_coordinates.x() );
  const auto z = world_coordinates.z();

  for( size_t itile = 0; itile < m_tiles.size(); ++itile )
  {
    const auto& tile = m_tiles.at(itile);

    if( std::abs( z - tile.m_centerZ ) > tile.m_sizeZ/2 ) continue;
    if( std::abs( bind_angle( phi - tile.m_centerPhi ) ) > tile.m_sizePhi/2 ) continue;

    return itile;
  }

  return -1;
}

//________________________________________________________________________________
int CylinderGeomMicromegas::find_strip_from_world_coords( uint tileid, const TVector3& world_coordinates ) const
{
  // convert to local coordinates
  const auto local_coordinates = get_local_from_world_coords( tileid, world_coordinates );
  return find_strip_from_local_coords( tileid, local_coordinates );
}

//________________________________________________________________________________
int CylinderGeomMicromegas::find_strip_from_local_coords( uint tileid, const TVector3& local_coordinates ) const
{

  // check against thickness
  if( std::abs( local_coordinates.y() ) > m_thickness/2 ) return -1;

  // get tile
  const auto& tile = m_tiles.at(tileid);

  // check azimuth
  if( std::abs( local_coordinates.x() ) > tile.m_sizePhi*reference_radius/2 ) return -1;
  
  // check z extend
  if( std::abs( local_coordinates.z() ) > tile.m_sizeZ/2 ) return -1;

  // we found a tile to which the hit belong
  // calculate strip index, depending on cylinder direction
  switch( m_segmentation_type )
  {
    case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
    return (int) std::floor( (local_coordinates.x() + tile.m_sizePhi*reference_radius/2)/m_pitch );

    case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
    return (int) std::floor( (local_coordinates.z() + tile.m_sizeZ/2)/m_pitch );
  }

  // unreachable
  return -1;
}

//________________________________________________________________________________
double CylinderGeomMicromegas::get_strip_length( uint tileid ) const
{
  assert( tileid < m_tiles.size() );
  switch( m_segmentation_type )
  {
    case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
    return m_tiles[tileid].m_sizeZ;

    case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
    return m_tiles[tileid].m_sizePhi*reference_radius;
  }

  // unreachable
  return 0;
}

//________________________________________________________________________________
uint CylinderGeomMicromegas::get_strip_count( uint tileid ) const
{
  assert( tileid < m_tiles.size() );
  switch( m_segmentation_type )
  {
    case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
    return std::floor( m_tiles[tileid].m_sizePhi*reference_radius/m_pitch );

    case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
    return std::floor( m_tiles[tileid].m_sizeZ/m_pitch );
  }

  // unreachable
  return 0;
}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_local_coordinates( uint tileid, uint stripnum ) const
{
  assert( tileid < m_tiles.size() );
  
  // get tile
  const auto& tile = m_tiles[tileid];
  
  switch( m_segmentation_type )
  {
    case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
    return TVector3( (0.5+stripnum)*m_pitch - tile.m_sizePhi*reference_radius/2, 0, 0 );
    
    case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
    return TVector3( 0, 0, (0.5+stripnum)*m_pitch - tile.m_sizeZ/2 );
  }
  
  // unreachable
  return TVector3();
}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_world_coordinates( uint tileid, uint stripnum ) const
{
  assert( tileid < m_tiles.size() );
  return get_world_from_local_coords( tileid, get_local_coordinates( tileid, stripnum ) );
}
 
//________________________________________________________________________________
void CylinderGeomMicromegas::identify( std::ostream& out ) const
{
  out << "CylinderGeomMicromegas" << std::endl;
  out << "layer: " << m_layer << std::endl;
  out << "segmentation_type: " << (m_segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_PHI ? "SEGMENTATION_PHI":"SEGMENTATION_Z") << std::endl;
  out << "drift_direction: " << (m_drift_direction == MicromegasDefs::DriftDirection::INWARD ? "INWARD":"OUTWARD") << std::endl;
  out << "radius: " << m_radius << "cm" << std::endl;
  out << "thickness: " << m_thickness << "cm" << std::endl;
  out << "zmin: " << m_zmin << "cm" << std::endl;
  out << "zmax: " << m_zmax << "cm" << std::endl;
  out << "pitch: " << m_pitch << "cm" << std::endl;
  out << std::endl;
}

//________________________________________________________________________________
bool  CylinderGeomMicromegas::check_radius( const TVector3& world_coordinates ) const
{
  const auto radius = get_r( world_coordinates.x(), world_coordinates.y() );
  return std::abs( radius - m_radius ) <= m_thickness/2;
}

//________________________________________________________________________________
TGeoHMatrix CylinderGeomMicromegas::transformation_matrix( uint tileid ) const
{
  assert( tileid < m_tiles.size() );
  const auto& tile = m_tiles[tileid];

  // local to master transformation matrix
  TGeoHMatrix matrix;

  // rotate around z axis
  matrix.RotateZ( tile.m_centerPhi*180./M_PI - 90 );

  // translate to tile center
  matrix.SetDx( m_radius*std::cos( tile.m_centerPhi ) );
  matrix.SetDy( m_radius*std::sin( tile.m_centerPhi ) );
  matrix.SetDz( tile.m_centerZ );

  return matrix;
}
