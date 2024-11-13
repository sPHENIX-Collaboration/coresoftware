#ifndef TPC_TPCGLOBALPOSITIONWRAPPER_H
#define TPC_TPCGLOBALPOSITIONWRAPPER_H

/*
 * \file TpcGlobalPositionWrapper.h
 * \brief provides the tpc 3d global position with all distortion and crossing corrections applied
 * \author Joe Osborn <josborn1@bnl.gov>, Hugo Pereira Da Costa <hugo.pereira-da-costa@lanl.gov>
 */
#include "TpcClusterZCrossingCorrection.h"
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

  //! load relevant nodes from tree
  void loadNodes(PHCompositeNode* /*topnode*/);

  //! apply all loaded distortion corrections to a given position
  Acts::Vector3 applyDistortionCorrections( Acts::Vector3 /*source*/ ) const;

  //! get distortion corrected global position from cluster
  /**
   * first converts cluster position local coordinate to global coordinates
   * then, for TPC clusters only, applies crossing correction, and distortion corrections
   */
  Acts::Vector3 getGlobalPositionDistortionCorrected(const TrkrDefs::cluskey&, TrkrCluster*, short int /*crossing*/ ) const;

  private:

  //! cluster z crossing correction interface
  TpcClusterZCrossingCorrection m_crossingCorrection;

  //! distortion correction interface
  TpcDistortionCorrection m_distortionCorrection;

  //! acts geometry
  ActsGeometry* m_tGeometry = nullptr;

  //! module edge distortion correction container
  TpcDistortionCorrectionContainer* m_dcc_module_edge{nullptr};

  //! static distortion correction container
  TpcDistortionCorrectionContainer* m_dcc_static{nullptr};

  //! average distortion correction container
  TpcDistortionCorrectionContainer* m_dcc_average{nullptr};

  //! fluctuation distortion container
  TpcDistortionCorrectionContainer* m_dcc_fluctuation{nullptr};

};

#endif
