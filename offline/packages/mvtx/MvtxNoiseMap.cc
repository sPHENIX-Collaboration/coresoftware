/**
 * @file mvtx/MvtxHitPruner.h
 * @author Yasser Corrales Morales
 * @date April 2024
 * @brief  Definition of the MVTX NoiseMap
 */

#include "MvtxNoiseMap.h"

ClassImp(MvtxNoiseMap);

//void NoiseMap::print()
//{
//  int nc = 0, np = 0, nm = 0;
//  for (const auto& map : mNoisyPixels) {
//    if (!map.empty()) {
//      nc++;
//    }
//    np += map.size();
//    if (map.find(getKey(-1, -1)) != map.end()) {
//      nm++;
//      nc--;
//      np--;
//    }
//  }
//  LOG(info) << "Number of fully maske chips " << nm;
//  LOG(info) << "Number of noisy chips: " << nc;
//  LOG(info) << "Number of noisy pixels: " << np;
//  LOG(info) << "Number of of strobes: " << mNumOfStrobes;
//  LOG(info) << "Probability threshold: " << mProbThreshold;
//}
//
//void NoiseMap::fill(const gsl::span<const CompClusterExt> data)
//{
//  for (const auto& c : data) {
//    if (c.getPatternID() != o2::itsmft::CompCluster::InvalidPatternID) {
//      // For the noise calibration, we use "pass1" clusters...
//      continue;
//    }
//
//    auto id = c.getSensorID();
//    auto row = c.getRow();
//    auto col = c.getCol();
//
//    // A simplified 1-pixel calibration
//    increaseNoiseCount(id, row, col);
//  }
//}
