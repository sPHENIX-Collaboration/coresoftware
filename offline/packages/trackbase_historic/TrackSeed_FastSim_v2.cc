#include "TrackSeed_FastSim_v2.h"

#include <iostream>
#include <map>

TrackSeed_FastSim_v2::TrackSeed_FastSim_v2(const TrackSeed& source)
{
  TrackSeed_FastSim_v2::CopyFrom(source);
}

void TrackSeed_FastSim_v2::CopyFrom(const TrackSeed& source)
{
  // do nothing if copying onto oneself
  if (this == &source)
  {
    return;
  }

  // parent class method
  TrackSeed_v2::CopyFrom(source);

  // additional members
  m_truth_track_id = source.get_truth_track_id();
  m_nmeas = source.get_num_measurements();
  m_g4hit_ids = source.g4hit_ids();
}

void TrackSeed_FastSim_v2::identify(std::ostream& os) const
{
  TrackSeed_v2::identify(os);

  os << "TrackSeed_FastSim_v2 Object ";
  os << "m_truth_Track_id: " << m_truth_track_id << std::endl;
  os << "m_nmeas: " << m_nmeas << std::endl;

  os << "G4Hit IDs:" << std::endl;
  for (const auto& pair : m_g4hit_ids)
  {
    os << "\thit container ID" << pair.first << " with hits: ";
    for (const auto& hitid : pair.second)
    {
      os << hitid << " ";
    }
    os << std::endl;
  }
  return;
}

float TrackSeed_FastSim_v2::get_phi(TrkrClusterContainer* /*clusters*/,
                                    ActsGeometry* /*tGeometry*/) const
{
  const auto [x, y] = findRoot();
  return std::atan2(-1 * (TrackSeed_v2::get_X0() - x), (TrackSeed_v2::get_Y0() - y));
}
