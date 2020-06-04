#ifndef MICROMEGASTILE_H
#define MICROMEGASTILE_H

/*!
 * \file MicromegasTile.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <vector>

//! header only class that contains information about a given Tile location inside CylinderGeom
class MicromegasTile
{

  public:

  using List = std::vector<MicromegasTile>;

  //! default constructor
  MicromegasTile()
  {}

  //! constructor
  MicromegasTile( double centerPhi, double centerZ, double sizePhi, double sizeZ )
    :m_centerPhi( centerPhi )
    ,m_centerZ( centerZ )
    ,m_sizePhi( sizePhi )
    ,m_sizeZ( sizeZ )
  {}

  double m_centerPhi = 0;
  double m_centerZ = 0;
  double m_sizePhi = 0;
  double m_sizeZ = 0;

};

#endif
