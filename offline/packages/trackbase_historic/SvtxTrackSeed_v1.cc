#include "SvtxTrackSeed_v1.h"


SvtxTrackSeed_v1::SvtxTrackSeed_v1()
{}

SvtxTrackSeed_v1::SvtxTrackSeed_v1(const SvtxTrackSeed& seed)
{ SvtxTrackSeed_v1::CopyFrom( seed ); }

SvtxTrackSeed_v1::SvtxTrackSeed_v1(const SvtxTrackSeed_v1& seed)
{ SvtxTrackSeed_v1::CopyFrom( seed ); }

SvtxTrackSeed_v1& SvtxTrackSeed_v1::operator=(const SvtxTrackSeed_v1& seed)
{ if( this != &seed ) CopyFrom( seed ); return *this; }

SvtxTrackSeed_v1::~SvtxTrackSeed_v1()
{}

void SvtxTrackSeed_v1::CopyFrom( const SvtxTrackSeed& seed )
{
  if( this == &seed ) return;
  SvtxTrackSeed::CopyFrom( seed );

  m_is_positive_charge = seed.get_positive_charge();
  m_px = seed.get_px();
  m_py = seed.get_py();
  m_pz = seed.get_pz();
  m_pcax = seed.get_x();
  m_pcay = seed.get_y();
  m_pcaz = seed.get_z();
  
  m_cluster_keys.clear();
  std::copy(seed.begin_cluster_keys(), seed.end_cluster_keys(), 
	    std::inserter(m_cluster_keys, m_cluster_keys.begin() ) );

}


void SvtxTrackSeed_v1::identify(std::ostream& os) const
{
  os << "SvtxTrackSeed_v1 object ";
  os << "charge " << get_charge() << std::endl;
  os << "(px,py,pz) = (" << get_px() << ", " << get_py()
     << ", " << get_pz() << ")" << std::endl;
  os << "(x,y,z) = (" << get_x() << ", " << get_y() << ", " << get_z() 
     << ")" << std::endl;
  os << "list of cluster keys ";
  if(m_cluster_keys.size() > 0) 
    {
      for (SvtxTrackSeed::ConstClusterKeyIter iter = begin_cluster_keys();
	   iter != end_cluster_keys();
	   ++iter)
	{
	  TrkrDefs::cluskey cluster_key = *iter;
	  os << cluster_key << ", ";
	}
    }
  
  os << std::endl;
  return;
}
