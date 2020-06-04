#ifndef MICROMEGASTILE_H
#define MICROMEGASTILE_H

/*!
 * \file MicromegasTile.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <array>
#include <vector>

//! header only class that contains information about a given Tile location inside CylinderGeom
class MicromegasTile
{

  public:

  using List = std::vector<MicromegasTile>;

  //! default constructor
  MicromegasTile()
  {}

  //! convenient construct from std::array
  MicromegasTile( const std::array<double, 4> values )
    :m_centerPhi( values[0] )
    ,m_centerZ( values[1] )
    ,m_sizePhi( values[2] )
    ,m_sizeZ( values[3] )
  {}

  double m_centerPhi = 0;
  double m_centerZ = 0;
  double m_sizePhi = 0;
  double m_sizeZ = 0;

};

#endif
