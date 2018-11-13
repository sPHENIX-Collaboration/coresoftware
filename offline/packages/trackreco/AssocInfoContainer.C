#include "../trackreco/AssocInfoContainer.h"

using namespace std;

ClassImp(AssocInfoContainer);

AssocInfoContainer::AssocInfoContainer() :
		_map_cluster_id_track_id()
{}

AssocInfoContainer::~AssocInfoContainer() {Reset();}

void AssocInfoContainer::Reset() {
	_map_cluster_id_track_id.clear();
}

void AssocInfoContainer::identify(std::ostream& os) const {
  cout << "---ClusterTrackMap--------------------------" << endl;
  for (auto iter = _map_cluster_id_track_id.begin();
  		iter != _map_cluster_id_track_id.end();
  		++iter) {
    cout
		<< "{" <<  iter->first
		<< " -> " << iter->second
		<< "}"
		<< endl;
  }
}
