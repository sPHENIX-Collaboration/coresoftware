#ifndef TPC_TPCDISTORTIONCORRECTION_H
#define TPC_TPCDISTORTIONCORRECTION_H

/*!
 * \file TpcDistortionCorrection.h
 * \brief applies provided distortion corrections to a cluster and returns corrected position
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <TVector3.h>

class TpcDistortionCorrectionObject;

class TpcDistortionCorrection
{
  public:
  
  //! constructor
  TpcDistortionCorrection() = default;
  
  //! get cluster corrected 3D position using given DistortionCorrectionObject
  TVector3 get_corrected_position( const TVector3&, const TpcDistortionCorrectionObject* ) const;
 
  private: 
  
  //! verbosity
  int m_verbosity = 0;
  
};

#endif
