#ifndef TRACKBASEHISTORIC_TRACKSEEDHELPER_H
#define TRACKBASEHISTORIC_TRACKSEEDHELPER_H

/*!
 * \file TrackSeedHelper.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@lanl.gov>
 */

#include <trackbase/TrkrDefs.h>

#include <Acts/Definitions/Algebra.hpp>

#include <cstdint>
#include <map>

class TrackSeed;

namespace TrackSeedHelper
{
  using position_map_t = std::map<TrkrDefs::cluskey, Acts::Vector3>;

  float get_phi(TrackSeed const*, const position_map_t&);
  float get_phi_fastsim(TrackSeed const*);
  void circleFitByTaubin(
    TrackSeed*, const position_map_t& positions,
    uint8_t startLayer = 0,
    uint8_t endLayer = 58);
  std::pair<float, float> findRoot(TrackSeed const* seed);
  std::pair<float, float> findRoot(const float& qOverR, const float& X0, const float& Y0);

  void lineFit(
      TrackSeed*, const position_map_t& positions,
      uint8_t startLayer = 0,
      uint8_t endLayer = 58);

  float get_x(TrackSeed const*);
  float get_y(TrackSeed const*);
  float get_z(TrackSeed const*);
  Acts::Vector3 get_xyz(TrackSeed const*);
}

#endif
