#ifndef WEIGHTED_TRACK_MAP_H
#define WEIGHTED_TRACK_MAP_H

#include "WeightedTrack.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <utility>

class WeightedTrackMap : public PHObject, private std::map<unsigned int, WeightedTrack*> {
public:
	WeightedTrackMap() = default;
	WeightedTrackMap(WeightedTrackMap const& /*other*/);
	WeightedTrackMap& operator=(WeightedTrackMap const& /*other*/);
 	~WeightedTrackMap() override { Reset(); }

	PHObject* CloneMe() const override { return new WeightedTrackMap(*this); }

 	void identify (std::ostream& /*out*/ = std::cout) const override;
	void Reset() override;
	int isValid() const override { return 1; }


	// Types
	using std::map<unsigned int, WeightedTrack*>::key_type;
	using std::map<unsigned int, WeightedTrack*>::mapped_type;
	using std::map<unsigned int, WeightedTrack*>::value_type;
	using std::map<unsigned int, WeightedTrack*>::size_type;
	using std::map<unsigned int, WeightedTrack*>::difference_type;
	using std::map<unsigned int, WeightedTrack*>::key_compare;
	using std::map<unsigned int, WeightedTrack*>::allocator_type;
	using std::map<unsigned int, WeightedTrack*>::reference;
	using std::map<unsigned int, WeightedTrack*>::const_reference;
	using std::map<unsigned int, WeightedTrack*>::pointer;
	using std::map<unsigned int, WeightedTrack*>::const_pointer;
	using std::map<unsigned int, WeightedTrack*>::iterator;
	using std::map<unsigned int, WeightedTrack*>::const_iterator;
	using std::map<unsigned int, WeightedTrack*>::reverse_iterator;
	using std::map<unsigned int, WeightedTrack*>::const_reverse_iterator;
	using std::map<unsigned int, WeightedTrack*>::node_type;
	using std::map<unsigned int, WeightedTrack*>::insert_return_type;

	// Element access
	using std::map<unsigned int, WeightedTrack*>::at;
	using std::map<unsigned int, WeightedTrack*>::operator[];

	// Iterators
	using std::map<unsigned int, WeightedTrack*>::begin;
	using std::map<unsigned int, WeightedTrack*>::cbegin;
	using std::map<unsigned int, WeightedTrack*>::end;
	using std::map<unsigned int, WeightedTrack*>::cend;
	using std::map<unsigned int, WeightedTrack*>::rbegin;
	using std::map<unsigned int, WeightedTrack*>::crbegin;
	using std::map<unsigned int, WeightedTrack*>::rend;
	using std::map<unsigned int, WeightedTrack*>::crend;

	// Capacity
	using std::map<unsigned int, WeightedTrack*>::empty;
	using std::map<unsigned int, WeightedTrack*>::size;
	using std::map<unsigned int, WeightedTrack*>::max_size;

	// Modifiers
	using std::map<unsigned int, WeightedTrack*>::clear;
	using std::map<unsigned int, WeightedTrack*>::insert;
	// using std::map<unsigned int, WeightedTrack*>::insert_range; c++23
	using std::map<unsigned int, WeightedTrack*>::insert_or_assign;
	using std::map<unsigned int, WeightedTrack*>::emplace;
	using std::map<unsigned int, WeightedTrack*>::emplace_hint;
	using std::map<unsigned int, WeightedTrack*>::try_emplace;
	using std::map<unsigned int, WeightedTrack*>::erase;
	using std::map<unsigned int, WeightedTrack*>::swap;
	using std::map<unsigned int, WeightedTrack*>::extract;
	using std::map<unsigned int, WeightedTrack*>::merge;

	// Lookup
	using std::map<unsigned int, WeightedTrack*>::count;
	using std::map<unsigned int, WeightedTrack*>::find;
	using std::map<unsigned int, WeightedTrack*>::contains;
	using std::map<unsigned int, WeightedTrack*>::equal_range;
	using std::map<unsigned int, WeightedTrack*>::lower_bound;
	using std::map<unsigned int, WeightedTrack*>::upper_bound;

private:
	ClassDefOverride(WeightedTrackMap, 1);
};

#endif//WEIGHTED_TRACK_MAP_H
