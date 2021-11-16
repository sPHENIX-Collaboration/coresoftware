#ifndef TPC_TPCDISTORTIONCORRECTIONOBJECT_H
#define TPC_TPCDISTORTIONCORRECTIONOBJECT_H

/*!
 * \file TpcDistortionCorrectionObject.h
 * \brief stores distortion correction histograms on the node tree
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <array>

class TH3;

class TpcDistortionCorrectionObject
{
  public:

  //! constructor
  TpcDistortionCorrectionObject() = default;
  
  //!@name space charge distortion histograms
  //@{
  std::array<TH3*,2> m_hDRint = {{nullptr, nullptr}};
  std::array<TH3*,2> m_hDPint = {{nullptr, nullptr}};
  std::array<TH3*,2> m_hDZint = {{nullptr, nullptr}};
  //@}

};

#endif
