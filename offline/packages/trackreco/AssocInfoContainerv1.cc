#include "AssocInfoContainerv1.h"

AssocInfoContainerv1::AssocInfoContainerv1()
  : _map_cluster_id_track_id()
{
}

AssocInfoContainerv1::~AssocInfoContainerv1() { Reset(); }

void AssocInfoContainerv1::Reset()
{
  _map_cluster_id_track_id.clear();
}

void AssocInfoContainerv1::identify(std::ostream& os) const
{
  os << "---ClusterTrackMap--------------------------" << std::endl;
  for (auto iter = _map_cluster_id_track_id.begin();
       iter != _map_cluster_id_track_id.end();
       ++iter)
  {
    os
        << "{" << iter->first
        << " -> " << iter->second
        << "}"
        << std::endl;
  }
}
