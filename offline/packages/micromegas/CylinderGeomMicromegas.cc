/*!
 * \file CylinderGeomMicromegas.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "CylinderGeomMicromegas.h"

#include <g4main/PHG4Hit.h>
#include <trackbase/ActsGeometry.h>

#include <Acts/Surfaces/RectangleBounds.hpp>
#include <Acts/Surfaces/SurfaceBounds.hpp>
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
TVector3 CylinderGeomMicromegas::get_local_from_world_coords( uint tileid, ActsGeometry* geometry, const TVector3& world_coordinates ) const
{
  assert( tileid < m_tiles.size() );
  const auto hitsetkey = MicromegasDefs::genHitSetKey(get_layer(), get_segmentation_type(), tileid);
  const auto surface = geometry->maps().getMMSurface(hitsetkey);

  const Acts::Vector3 global(
    world_coordinates.x()*Acts::UnitConstants::cm,
    world_coordinates.y()*Acts::UnitConstants::cm,
    world_coordinates.z()*Acts::UnitConstants::cm );

  // convert to local
  /* this is equivalent to calling surface->globalToLocal but without the "on surface" check, and while returning a full Acts::Vector3 */
  const auto local = surface->transform(geometry->geometry().getGeoContext()).inverse()*global;
  return TVector3(
    local.x()/Acts::UnitConstants::cm,
    local.y()/Acts::UnitConstants::cm,
    local.z()/Acts::UnitConstants::cm );
}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_local_from_world_vect( uint tileid, ActsGeometry* geometry, const TVector3& world_vect ) const
{
  assert( tileid < m_tiles.size() );
  const auto hitsetkey = MicromegasDefs::genHitSetKey(get_layer(), get_segmentation_type(), tileid);
  const auto surface = geometry->maps().getMMSurface(hitsetkey);

  const Acts::Vector3 global(
    world_vect.x()*Acts::UnitConstants::cm,
    world_vect.y()*Acts::UnitConstants::cm,
    world_vect.z()*Acts::UnitConstants::cm );

  const Acts::Vector3 local = surface->referenceFrame(geometry->geometry().getGeoContext(), Acts::Vector3(), Acts::Vector3()).inverse()*global;
  return TVector3(
    local.x()/Acts::UnitConstants::cm,
    local.y()/Acts::UnitConstants::cm,
    local.z()/Acts::UnitConstants::cm );
}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_world_from_local_coords( uint tileid, ActsGeometry* geometry, const TVector2& local_coordinates ) const
{
  assert( tileid < m_tiles.size() );
  const auto hitsetkey = MicromegasDefs::genHitSetKey(get_layer(), get_segmentation_type(), tileid);
  const auto surface = geometry->maps().getMMSurface(hitsetkey);

  const Acts::Vector2 local(
    local_coordinates.X()*Acts::UnitConstants::cm,
    local_coordinates.Y()*Acts::UnitConstants::cm );

  // convert to global
  const auto global =  surface->localToGlobal(geometry->geometry().getGeoContext(), local, Acts::Vector3());
  return TVector3(
    global.x()/Acts::UnitConstants::cm,
    global.y()/Acts::UnitConstants::cm,
    global.z()/Acts::UnitConstants::cm );
}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_world_from_local_coords( uint tileid, ActsGeometry* geometry, const TVector3& local_coordinates ) const
{
  assert( tileid < m_tiles.size() );
  const auto hitsetkey = MicromegasDefs::genHitSetKey(get_layer(), get_segmentation_type(), tileid);
  const auto surface = geometry->maps().getMMSurface(hitsetkey);

  const Acts::Vector3 local(
    local_coordinates.x()*Acts::UnitConstants::cm,
    local_coordinates.y()*Acts::UnitConstants::cm,
    local_coordinates.z()*Acts::UnitConstants::cm );

  // convert to global
  /* this is equivalent to calling surface->localToGlobal but without assuming that the local point is on surface */
  const auto global =  surface->transform(geometry->geometry().getGeoContext())*local;
  return TVector3(
    global.x()/Acts::UnitConstants::cm,
    global.y()/Acts::UnitConstants::cm,
    global.z()/Acts::UnitConstants::cm );
}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_world_from_local_vect( uint tileid, ActsGeometry* geometry, const TVector3& local_vect ) const
{
  assert( tileid < m_tiles.size() );
  const auto hitsetkey = MicromegasDefs::genHitSetKey(get_layer(), get_segmentation_type(), tileid);
  const auto surface = geometry->maps().getMMSurface(hitsetkey);

  Acts::Vector3 local(
    local_vect.x()*Acts::UnitConstants::cm,
    local_vect.y()*Acts::UnitConstants::cm,
    local_vect.z()*Acts::UnitConstants::cm );

  const Acts::Vector3 global = surface->referenceFrame(geometry->geometry().getGeoContext(), Acts::Vector3(), Acts::Vector3())*local;
  return TVector3(
    global.x()/Acts::UnitConstants::cm,
    global.y()/Acts::UnitConstants::cm,
    global.z()/Acts::UnitConstants::cm );
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
int CylinderGeomMicromegas::find_strip_from_world_coords( uint tileid, ActsGeometry* geometry, const TVector3& world_coordinates ) const
{
  // convert to local coordinates
  const auto local_coordinates = get_local_from_world_coords( tileid, geometry, world_coordinates );
  
  // check against thickness
  if( std::abs(local_coordinates.z() ) >= m_thickness/2 )
  {
    std::cout << " CylinderGeomMicromegas::find_strip_from_world_coords - invalid world coordinates" << std::endl;
    return -1;
  } else {
    return find_strip_from_local_coords( tileid, geometry, { local_coordinates.x(), local_coordinates.y() } );
  }
}

//________________________________________________________________________________
int CylinderGeomMicromegas::find_strip_from_local_coords( uint tileid, ActsGeometry* geometry, const TVector2& local_coordinates ) const
{
  assert( tileid < m_tiles.size() );
  const auto hitsetkey = MicromegasDefs::genHitSetKey(get_layer(), get_segmentation_type(), tileid);
  const auto surface = geometry->maps().getMMSurface(hitsetkey);

  // get boundaries and corresponding dimension
  assert( surface->bounds().type() == Acts::SurfaceBounds::BoundsType::eRectangle );  
  auto rectangle_bounds = static_cast<const Acts::RectangleBounds*>( &surface->bounds() );
  const auto half_length_x = rectangle_bounds->halfLengthX()/Acts::UnitConstants::cm;
  const auto half_length_y = rectangle_bounds->halfLengthY()/Acts::UnitConstants::cm;

  // check azimuth
  if( std::abs( local_coordinates.X() ) >  half_length_x ) return -1;

  // check z extend
  if( std::abs( local_coordinates.Y() ) > half_length_y ) return -1;

  // calculate strip index, depending on segmentation
  switch( m_segmentation_type )
  {
    case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
    return (int) std::floor((local_coordinates.X() + half_length_x)/m_pitch);

    case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
    return (int) std::floor((local_coordinates.Y() + half_length_y)/m_pitch);
  }

  // unreachable
  return -1;
}

//________________________________________________________________________________
double CylinderGeomMicromegas::get_strip_length( uint tileid, ActsGeometry* geometry ) const
{
  assert( tileid < m_tiles.size() );
  const auto hitsetkey = MicromegasDefs::genHitSetKey(get_layer(), get_segmentation_type(), tileid);
  const auto surface = geometry->maps().getMMSurface(hitsetkey);

  // get boundaries and return corresponding dimension
  assert( surface->bounds().type() == Acts::SurfaceBounds::BoundsType::eRectangle );  
  auto rectangle_bounds = static_cast<const Acts::RectangleBounds*>( &surface->bounds() );
  switch( m_segmentation_type )
  {
    case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
    return 2.*rectangle_bounds->halfLengthY()/Acts::UnitConstants::cm;

    case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
    return 2.*rectangle_bounds->halfLengthX()/Acts::UnitConstants::cm;
  }

  // unreachable
  return 0;
}

//________________________________________________________________________________
uint CylinderGeomMicromegas::get_strip_count( uint tileid, ActsGeometry* geometry ) const
{
  assert( tileid < m_tiles.size() );
  const auto hitsetkey = MicromegasDefs::genHitSetKey(get_layer(), get_segmentation_type(), tileid);
  const auto surface = geometry->maps().getMMSurface(hitsetkey);

  // get boundaries and corresponding dimension
  assert( surface->bounds().type() == Acts::SurfaceBounds::BoundsType::eRectangle );  
  auto rectangle_bounds = static_cast<const Acts::RectangleBounds*>( &surface->bounds() );
  switch( m_segmentation_type )
  {
    case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
    return std::floor( (2.*rectangle_bounds->halfLengthX()/Acts::UnitConstants::cm)/m_pitch );

    case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
    return std::floor( (2.*rectangle_bounds->halfLengthY()/Acts::UnitConstants::cm)/m_pitch );
  }

  // unreachable
  return 0;
}

//________________________________________________________________________________
TVector2 CylinderGeomMicromegas::get_local_coordinates( uint tileid, ActsGeometry* geometry, uint stripnum ) const
{
  assert( tileid < m_tiles.size() );
  const auto hitsetkey = MicromegasDefs::genHitSetKey(get_layer(), get_segmentation_type(), tileid);
  const auto surface = geometry->maps().getMMSurface(hitsetkey);

  // get boundaries and return corresponding dimension
  assert( surface->bounds().type() == Acts::SurfaceBounds::BoundsType::eRectangle );  
  auto rectangle_bounds = static_cast<const Acts::RectangleBounds*>( &surface->bounds() );

  switch( m_segmentation_type )
  {
    case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
    return TVector2( (0.5+stripnum)*m_pitch-rectangle_bounds->halfLengthX()/Acts::UnitConstants::cm, 0 );

    case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
    return TVector2( 0, (0.5+stripnum)*m_pitch-rectangle_bounds->halfLengthY()/Acts::UnitConstants::cm );
  }

  // unreachable
  return TVector2();
}

//________________________________________________________________________________
TVector3 CylinderGeomMicromegas::get_world_coordinates( uint tileid, ActsGeometry* geometry, uint stripnum ) const
{ return get_world_from_local_coords( tileid, geometry, get_local_coordinates( tileid, geometry, stripnum ) ); }

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
