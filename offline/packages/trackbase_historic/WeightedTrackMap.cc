#include "WeightedTrackMap.h"

WeightedTrackMap::WeightedTrackMap (
	WeightedTrackMap const& other
) : PHObject(other), std::map<unsigned int, WeightedTrack*>(other) {
	for (auto const& [track_id, weighted_track] : other) {
		insert({track_id, dynamic_cast<WeightedTrack*>(weighted_track->CloneMe())});
	}
}

WeightedTrackMap&
WeightedTrackMap::operator= (
	WeightedTrackMap const& other
) {
	Reset();
	for (auto const& [track_id, weighted_track] : other) {
		insert({track_id, dynamic_cast<WeightedTrack*>(weighted_track->CloneMe())});
	}
	return *this;
}

void
WeightedTrackMap::Reset (
) {
	for (auto& [track_id, weighted_track] : *this) {
		delete weighted_track;
	}
	clear();
}

void
WeightedTrackMap::identify (
	std::ostream& out
) const {
	out
		<< "WeightedTrackMap "
		<< "size: " << size()
		<< std::endl;
}
