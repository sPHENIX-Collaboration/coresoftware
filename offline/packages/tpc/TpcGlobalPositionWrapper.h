#ifndef TPC_TPCGLOBALPOSITIONWRAPPER_H
#define TPC_TPCGLOBALPOSITIONWRAPPER_H

/*
 * \file TpcGlobalPositionWrapper.h
 * \brief provides the tpc 3d global position with all distortion and crossing corrections applied
 * \author Joe Osborn <josborn1@bnl.gov>, Hugo Pereira Da Costa <hugo.pereira-da-costa@lanl.gov>
 */
#include "TpcDistortionCorrection.h"

#include <trackbase/TrkrDefs.h>


class ActsGeometry;
class PHCompositeNode;
class TpcDistortionCorrectionContainer;
class TrkrCluster;

class TpcGlobalPositionWrapper
{
  public:

  //! constructor
  explicit TpcGlobalPositionWrapper() = default;

  //! verbosity
  void set_verbosity(int value)
  {
    m_verbosity = value;
  }
  
  void set_suppressCrossing(bool value)
  {
    m_suppressCrossing = value;
  }

  //! verbosity
  int verbosity() const
  {
    return m_verbosity;
  }

  //! load relevant nodes from tree
  void loadNodes(PHCompositeNode* /*topnode*/);

  void set_enable_module_edge_corr(bool flag) { m_enable_module_edge_corr = flag; }
  void set_enable_static_corr(bool flag) { m_enable_static_corr = flag; }
  void set_enable_average_corr(bool flag) { m_enable_average_corr = flag; }
  void set_enable_fluctuation_corr(bool flag) { m_enable_fluctuation_corr = flag; }

  //! apply all loaded distortion corrections to a given position
  Acts::Vector3 applyDistortionCorrections( Acts::Vector3 /*source*/ ) const;

  //! get distortion corrected global position from cluster
  /**
   * first converts cluster position local coordinate to global coordinates
   * then, for TPC clusters only, applies crossing correction, and distortion corrections
   */
  Acts::Vector3 getGlobalPositionDistortionCorrected(const TrkrDefs::cluskey&, TrkrCluster*, short int /*crossing*/ ) const;

  private:

  //! verbosity
  unsigned int m_verbosity = 0;

  bool m_suppressCrossing = false;

  //! distortion correction interface
  TpcDistortionCorrection m_distortionCorrection;

  //! acts geometry
  ActsGeometry* m_tGeometry = nullptr;

  //! module edge distortion correction container
  TpcDistortionCorrectionContainer* m_dcc_module_edge{nullptr};
  bool m_enable_module_edge_corr = true;

  //! static distortion correction container
  TpcDistortionCorrectionContainer* m_dcc_static{nullptr};
  bool m_enable_static_corr = true;

  //! average distortion correction container
  TpcDistortionCorrectionContainer* m_dcc_average{nullptr};
  bool m_enable_average_corr = true;

  //! fluctuation distortion container
  TpcDistortionCorrectionContainer* m_dcc_fluctuation{nullptr};
  bool m_enable_fluctuation_corr = true;

};

#endif
