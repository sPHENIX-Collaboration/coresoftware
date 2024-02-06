#include "CentralityInfov2.h"

#include <limits>

void CentralityInfov2::identify(std::ostream& os) const
{
  os << "  CentralityInfo: " << std::endl;
  os << "      Centrile: " << (has_centile(CentralityInfo::PROP::mbd_NS) ? get_centile(CentralityInfo::PROP::mbd_NS) : -999.99) << std::endl;
  os << "      CentBin : " << (has_centrality_bin(CentralityInfo::PROP::mbd_NS) ? get_centrality_bin(CentralityInfo::PROP::mbd_NS) : -999.99) << std::endl;

  return;
}

bool CentralityInfov2::has_centrality_bin(const PROP prop_id) const
{
  return _centrality_bin_map.find(prop_id) != _centrality_bin_map.end();
}

void CentralityInfov2::set_centrality_bin(const PROP prop_id, int value)
{
  _centrality_bin_map[prop_id] = value;
}

int CentralityInfov2::get_centrality_bin(const PROP prop_id) const
{
  if (!has_centrality_bin(prop_id))
  {
    return std::numeric_limits<int>::quiet_NaN();
  }
  else
  {
    return _centrality_bin_map.at(prop_id);
  }
}
