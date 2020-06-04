/*!
 * \file CylinderGeomMicromegas.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "CylinderGeomMicromegas.h"

#include <TVector3.h>

namespace
{
  // convenient square function
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  // bind angle to [-M_PI,+M_PI[. This is useful when making the difference between two azimuthal angles
  template<class T>
    inline constexpr T bind_angle( const T& angle )
  {
    if( angle >= M_PI ) return angle - 2*M_PI;
    else if( angle < -M_PI ) return angle + 2*M_PI;
    else return angle;
  }

}

//________________________________________________________________________________
void CylinderGeomMicromegas::identify( std::ostream& out ) const
{
  out << "CylinderGeomMicromegas" << std::endl;
  out << "layer: " << m_layer << std::endl;
  out << "segmentation_type: " << (m_segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_PHI ? "SEGMENTATION_PHI":"SEGMENTATION_Z") << std::endl;
  out << "radius: " << m_radius << "cm" << std::endl;
  out << "thickness: " << m_thickness << "cm" << std::endl;
  out << "zmin: " << m_zmin << "cm" << std::endl;
  out << "zmax: " << m_zmax << "cm" << std::endl;
  out << std::endl;
}

//________________________________________________________________________________
CylinderGeomMicromegas::StripId CylinderGeomMicromegas::find_strip( const TVector3& world_location ) const
{

  // convert to polar coordinates
  const auto radius = std::sqrt( square( world_location.x() ) + square( world_location.y() ) );
  const auto phi = std::atan2( world_location.y(), world_location.x() );
  const auto z = world_location.z();

  if( std::abs( radius - m_radius ) > m_thickness/2 ) return std::make_pair( -1, -1 );

  for( size_t itile = 0; itile < m_tiles.size(); ++itile )
  {
    const auto& tile = m_tiles.at(itile);

    if( std::abs( z - tile.m_centerZ ) > tile.m_sizeZ/2 ) continue;
    if( std::abs( bind_angle( phi - tile.m_centerPhi ) ) > tile.m_sizePhi/2 ) continue;

    // we found a tile to which the hit belong
    // calculate strip index, depending on cylinder direction
    switch( m_segmentation_type )
    {

      case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
      return std::make_pair( itile, (int) std::floor( (bind_angle( phi - tile.m_centerPhi ) + tile.m_sizePhi/2)/tile.m_pitch ) );

      case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
      return std::make_pair( itile, (int) std::floor( (z - tile.m_centerZ + tile.m_sizeZ/2)/tile.m_pitch ) );

    }

  }

  // no tile found
  return  std::make_pair( -1, -1 );

}
