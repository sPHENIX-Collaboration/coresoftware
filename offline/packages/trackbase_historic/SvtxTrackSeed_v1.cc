#include "SvtxTrackSeed_v1.h"

SvtxTrackSeed_v1::SvtxTrackSeed_v1()
{}

// have to suppress uninitMenberVar from cppcheck since it triggers many false positive
// cppcheck-suppress uninitMemberVar
SvtxTrackSeed_v1::SvtxTrackSeed_v1(const SvtxTrackSeed_v1& seed)
{ SvtxTrackSeed_v1::CopyFrom( seed ); }

SvtxTrackSeed_v1& SvtxTrackSeed_v1::operator=(const SvtxTrackSeed_v1& seed)
{ if( this != &seed ) CopyFrom( seed ); return *this; }

SvtxTrackSeed_v1::~SvtxTrackSeed_v1()
{}

void SvtxTrackSeed_v1::CopyFrom( const TrackSeed& seed )
{
  if( this == &seed ) return;
  TrackSeed::CopyFrom( seed );

  m_silicon_seed = seed.get_silicon_seed_index();
  m_tpc_seed = seed.get_tpc_seed_index();

}

void SvtxTrackSeed_v1::identify(std::ostream& os) const
{

  os << "SvtxTrackSeed_v1 object ";
  os << "Silicon seed index: " << m_silicon_seed;
  os << "tpc seed index: " << m_tpc_seed;
}
