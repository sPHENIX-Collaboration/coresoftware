#include "ClusterHitsVerbose_v1.h"

std::vector<std::pair<int,float>>& ClusterHitsVerbose_v1::get_phi_hits()
{ return m_phi_hits; }

std::vector<std::pair<int,float>>& ClusterHitsVerbose_v1::get_z_hits()
{ return m_z_hits; }
