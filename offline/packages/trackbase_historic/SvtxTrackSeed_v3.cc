#include "SvtxTrackSeed_v3.h"

SvtxTrackSeed_v3::SvtxTrackSeed_v3() = default;

// have to suppress missingMemberCopy from cppcheck, it does not
// go down to the CopyFrom method where things are done correctly
// cppcheck-suppress missingMemberCopy
SvtxTrackSeed_v3::SvtxTrackSeed_v3(const SvtxTrackSeed_v3& seed)
  : TrackSeed(seed)
{
  SvtxTrackSeed_v3::CopyFrom(seed);
}

SvtxTrackSeed_v3& SvtxTrackSeed_v3::operator=(const SvtxTrackSeed_v3& seed)
{
  if (this != &seed)
  {
    CopyFrom(seed);
  }
  return *this;
}

SvtxTrackSeed_v3::~SvtxTrackSeed_v3() = default;

void SvtxTrackSeed_v3::CopyFrom(const TrackSeed& seed)
{
  if (this == &seed)
  {
    return;
  }
  TrackSeed::CopyFrom(seed);

  m_silicon_seed = seed.get_silicon_seed_index();
  m_tpc_seed = seed.get_tpc_seed_index();
  m_crossing_estimate = seed.get_crossing_estimate();
  m_crossing = seed.get_crossing();
}

void SvtxTrackSeed_v3::identify(std::ostream& os) const
{
  os << "SvtxTrackSeed_v3 object ";
  os << " Silicon seed index: " << m_silicon_seed;
  os << " tpc seed index: " << m_tpc_seed;
  os << " crossing estimate: " << m_crossing_estimate;
  os << " crossing: " << m_crossing;
  os << std::endl;
}
