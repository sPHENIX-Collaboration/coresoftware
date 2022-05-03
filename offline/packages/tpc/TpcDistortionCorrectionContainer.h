#ifndef TPC_TPCDISTORTIONCORRECTIONCONTAINER_H
#define TPC_TPCDISTORTIONCORRECTIONCONTAINER_H

/*!
 * \file TpcDistortionCorrectionContainer.h
 * \brief stores distortion correction histograms on the node tree
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <array>

class TH1;

class TpcDistortionCorrectionContainer
{
  public:

  //! constructor
  TpcDistortionCorrectionContainer() = default;

  //! flag to tell us whether to read z data or just interpolate
  int dimensions=3;
  
  //!@name space charge distortion histograms
  //@{
  std::array<TH1*,2> m_hDRint = {{nullptr, nullptr}};
  std::array<TH1*,2> m_hDPint = {{nullptr, nullptr}};
  std::array<TH1*,2> m_hDZint = {{nullptr, nullptr}};

  /// keep track of number of entries in each bin
  /** 
   * used temporarily  when building distortion corrections on the fly
   * it is not used to actually apply corrections to a given 3D point 
   */
  std::array<TH1*,2> m_hentries = {{nullptr, nullptr}};
  //@}


};

#endif
