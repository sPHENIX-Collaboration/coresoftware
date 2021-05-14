#include "AssocInfoContainerv1.h"

using namespace std;

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
  cout << "---ClusterTrackMap--------------------------" << endl;
  for (auto iter = _map_cluster_id_track_id.begin();
       iter != _map_cluster_id_track_id.end();
       ++iter)
  {
    cout
        << "{" << iter->first
        << " -> " << iter->second
        << "}"
        << endl;
  }
}
