// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef MICROMEGAS_MICROMEGASTILE_H
#define MICROMEGAS_MICROMEGASTILE_H

/*!
 * \file MicromegasTile.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <phool/PHObject.h>

#include <array>
#include <cassert>
#include <vector>

//! header only class that contains information about a given Tile location inside CylinderGeom
class MicromegasTile: public PHObject
{

  public:

  using List = std::vector<MicromegasTile>;

  //! default constructor
  MicromegasTile() = default;

  //! destructor
  ~MicromegasTile() override = default;

  //! constructor
  MicromegasTile( std::array<double, 4> values )
    :m_centerPhi( values[0] )
    ,m_centerZ( values[1] )
    ,m_sizePhi( values[2] )
    ,m_sizeZ( values[3] )
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

  ClassDefOverride(MicromegasTile,1)
};

#endif
