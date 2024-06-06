#ifndef TPC_TPCGLOBALPOSITIONWRAPPER_H
#define TPC_TPCGLOBALPOSITIONWRAPPER_H

/*
 * \file TpcGlobalPositionWrapper.h
 * \brief provides the tpc 3d global position with all distortion and crossing corrections applied
 * \author Joe Osborn <josborn1@bnl.gov>
 */
#include "TpcClusterZCrossingCorrection.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>
#include <Acts/Definitions/Algebra.hpp>

class TpcDistortionCorrection;
class TpcDistortionCorrectionContainer;
class TrkrCluster;

class TpcGlobalPositionWrapper
{
 public:
  static Acts::Vector3 getGlobalPositionDistortionCorrected(const TrkrDefs::cluskey& key, TrkrCluster* cluster, ActsGeometry* tGeometry, short int crossing, const TpcDistortionCorrectionContainer* staticCorrection, const TpcDistortionCorrectionContainer* averageCorrection, const TpcDistortionCorrectionContainer* fluctuationCorrection);
};

#endif
