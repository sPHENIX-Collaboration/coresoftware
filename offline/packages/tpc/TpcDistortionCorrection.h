#ifndef TPC_TPCDISTORTIONCORRECTION_H
#define TPC_TPCDISTORTIONCORRECTION_H

/*!
 * \file TpcDistortionCorrection.h
 * \brief applies provided distortion corrections to a 3D point and returns corrected position
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <Acts/Definitions/Algebra.hpp>

class TpcDistortionCorrectionContainer;

class TpcDistortionCorrection
{
  public:
  
  //! constructor
  TpcDistortionCorrection() = default;

  enum DistortionType
  {
    StaticOnly=0,
    BeamInducedAverage=1,
    BeamInducedFluctuation=2
  };
  
  enum CoordMask
  {
    COORD_PHI = 1<<0,
    COORD_R = 1<<1,
    COORD_Z = 1<<2,
    COORD_PHIZ = COORD_PHI|COORD_Z,
    COORD_ALL = COORD_PHI|COORD_R|COORD_Z
  };
  
  //! get cluster corrected 3D position using given DistortionCorrectionObject
 Acts::Vector3 get_corrected_position( const Acts::Vector3&, const TpcDistortionCorrectionContainer*, 
					unsigned int mask = COORD_ALL ) const;

};

#endif
