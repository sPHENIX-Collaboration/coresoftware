#include "CosmicTrackSeed_v1.h"

CosmicTrackSeed_v1::CosmicTrackSeed_v1()
{
}

// have to suppress uninitMenberVar from cppcheck since it triggers many false positive
// cppcheck-suppress uninitMemberVar
CosmicTrackSeed_v1::CosmicTrackSeed_v1(const CosmicTrackSeed_v1& seed)
{
  CosmicTrackSeed_v1::CopyFrom(seed);
}

CosmicTrackSeed_v1& CosmicTrackSeed_v1::operator=(const CosmicTrackSeed_v1& seed)
{
  if (this != &seed) CopyFrom(seed);
  return *this;
}

CosmicTrackSeed_v1::~CosmicTrackSeed_v1()
{
}

void CosmicTrackSeed_v1::CopyFrom(const TrackSeed& seed)
{
  if (this == &seed) return;
  TrackSeed::CopyFrom(seed);

  m_silicon_seed1 = seed.get_silicon_seed_index();
  m_tpc_seed1 = seed.get_tpc_seed_index();
  m_silicon_seed2 = seed.get_silicon_seed_index2();
  m_tpc_seed2 = seed.get_tpc_seed_index2();
}

void CosmicTrackSeed_v1::identify(std::ostream& os) const
{
  os << "CosmicTrackSeed_v1 object ";
  os << "Silicon seed index1: " << m_silicon_seed1;
  os << "tpc seed index1: " << m_tpc_seed1;
  os << "Silicon seed index2: " << m_silicon_seed2;
  os << "tpc seed index2: " << m_tpc_seed2;
}
