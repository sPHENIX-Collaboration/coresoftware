
#include "EvalLinksV1.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <set>
#include <map>
#include <float.h>
#include <limits.h>

ClassImp(EvalLinksV1)

using namespace std;

// constant used to record dangling links
// needed for forward evaluation to record
// missing reconstruction evals
// e.g. (g4particle,nothing) => nclusters
// or noise generated evals
// e.g. (cluster, noise) => nhits
const unsigned int EvalLinksV1::NULLID = UINT_MAX;

EvalLinksV1::EvalLinksV1(const std::string& left_name,
			 const std::string& right_name,
			 const std::string& weight_name) :
  EvalLinks(left_name,right_name,weight_name) {
  _stale = true;
  _left_name = left_name;
  _right_name = right_name;
  _weight_name = weight_name;
  _links.clear();
  _left_right_mmap.clear();
  _right_left_mmap.clear();
  _left_right_map.clear();
  _right_left_map.clear();
}

void EvalLinksV1::identify(std::ostream& os) const {
  if (stale()) refresh();

  os << "---EvalLinksV1--------------------------" << endl;
  os << " left:  " << _left_name << endl;
  os << " right: " << _right_name << endl;
  os << " weight: " << _weight_name << endl;

  os << " size = " << size() << endl;

  os << " links: " << endl;
  os << "   " << _left_name << " id <==> " << _right_name << " id : " << _weight_name << endl;
  os << "   " << "-------------------------------" << endl;
  for (std::map< std::pair<unsigned int,unsigned int>, float>::const_iterator citer = _links.begin();
       citer != _links.end();
       ++citer) {
    unsigned int left_id = citer->first.first;
    unsigned int right_id = citer->first.second;
    float weight = citer->second;
    
    os << "   " << left_id << " <==> " << right_id << " : " << weight;
    if (right_id == max_right(left_id)) os << " *";
    os << endl;
  }
  
  return;
}

void EvalLinksV1::set_names(const std::string& left_name,
			    const std::string& right_name,
			    const std::string& weight_name) {
  _left_name = left_name;
  _right_name = right_name;
  _weight_name = weight_name;
}

void EvalLinksV1::link(unsigned int left_id, 
		       unsigned int right_id,
		       float weight) {
  if (stale()) refresh();

  // insert into storage
  _links.insert(make_pair(make_pair(left_id,right_id),weight));

  // add link to multimaps
  _left_right_mmap.insert(make_pair(left_id,right_id));
  _right_left_mmap.insert(make_pair(right_id,left_id));

  // recalculate max links
  calc_max_right(left_id);
  calc_max_left(right_id);

  return;
}

void EvalLinksV1::unlink(unsigned int left_id, unsigned int right_id) {
  if (stale()) refresh();

  if (_links.find(make_pair(left_id,right_id)) != _links.end()) {
    _links.erase(make_pair(left_id,right_id));
    refresh();
  }
  
  return;
}

void EvalLinksV1::unlink_subleading() {
  if (stale()) refresh();

  for (std::map< std::pair<unsigned int,unsigned int>, float>::const_iterator citer = _links.begin();
       citer != _links.end();
       ++citer) {
    unsigned int left_id = citer->first.first;
    unsigned int right_id = citer->first.second;
    unsigned int max_right_id = max_right(left_id);

    if (right_id != max_right_id) {
      unlink(left_id,right_id);
    }    
  }
  
  return;
}

void EvalLinksV1::clear() {
  _stale = true;
  _left_name.clear();
  _right_name.clear();
  _weight_name.clear();
  _links.clear();
  _left_right_mmap.clear();
  _right_left_mmap.clear();
  _left_right_map.clear();
  _right_left_map.clear();
}


size_t EvalLinksV1::size() const {
  if (stale()) refresh();
  return _links.size();
}

bool EvalLinksV1::has_link(unsigned int left_id, unsigned int right_id) const {
  if (stale()) refresh();

  if (_links.find(make_pair(left_id,right_id)) != _links.end()) {
    return true;
  }

  return false;
}

float EvalLinksV1::get_weight(unsigned int left_id, unsigned int right_id) const {
  if (stale()) refresh();

  std::pair<unsigned int,unsigned int> pair = make_pair(left_id,right_id);
  std::map<std::pair<unsigned int,unsigned int>, float >::const_iterator citer = _links.find(pair);
  if (citer != _links.end()) {
    return citer->second;
  }

  return NAN;
}

std::set<unsigned int> EvalLinksV1::left(unsigned int right_id) const {
  if (stale()) refresh();
  
  std::set<unsigned int> left_links;
  for (std::multimap<unsigned int,unsigned int>::const_iterator citer = _right_left_mmap.lower_bound(right_id);
       citer != _right_left_mmap.upper_bound(right_id);
       ++citer) {
    left_links.insert(citer->second);
  }
  
  return left_links;
}

std::set<unsigned int> EvalLinksV1::right(unsigned int left_id) const {
  if (stale()) refresh();
   
  std::set<unsigned int> right_links;
  for (std::multimap<unsigned int,unsigned int>::const_iterator citer = _left_right_mmap.lower_bound(left_id);
       citer != _left_right_mmap.upper_bound(left_id);
       ++citer) {
    right_links.insert(citer->second);
  }
  
  return right_links;
}

unsigned int EvalLinksV1::max_left(unsigned int right_id) const {
  if (stale()) refresh();
  return _right_left_map[right_id];
}

unsigned int EvalLinksV1::max_right(unsigned int left_id) const {
  if (stale()) refresh();
  return _left_right_map[left_id];
}

void EvalLinksV1::refresh() const {

  _left_right_mmap.clear();
  _right_left_mmap.clear();
  _left_right_map.clear();
  _right_left_map.clear();

  // loop over all link pairs and recreate multimaps
  for (std::map<std::pair<unsigned int,unsigned int>,float>::const_iterator citer = _links.begin();
       citer != _links.end();
       ++citer) {
    unsigned int left_id = citer->first.first;
    unsigned int right_id = citer->first.second;

    _left_right_mmap.insert(make_pair(left_id,right_id));
    _right_left_mmap.insert(make_pair(right_id,left_id));    
  }

  // loop over all multimap keys and recompute max links
  for (std::map<unsigned int,unsigned int>::const_iterator citer = _left_right_mmap.begin();
       citer != _left_right_mmap.end();
       ++citer) {
    unsigned int left_id = citer->first;
    calc_max_right(left_id);
  }

  for (std::map<unsigned int,unsigned int>::const_iterator citer = _right_left_mmap.begin();
       citer != _right_left_mmap.end();
       ++citer) {
    unsigned int right_id = citer->first;
    calc_max_left(right_id);
  }
  
  _stale = false;
}

void EvalLinksV1::calc_max_left(unsigned int right_id) const {
  
  std::set<unsigned int> candidates = left(right_id);

  unsigned int max_left = 0xFFFFFFF;
  float max_weight = FLT_MIN;

  for (std::set<unsigned int>::const_iterator iter = candidates.begin();
       iter != candidates.end();
       ++iter) {
    unsigned int candidate = *iter;
    float candidate_weight = get_weight(candidate,right_id);
    
    if (candidate_weight > max_weight) {
      max_left = candidate;
      max_weight = candidate_weight;
    }
  }

  _right_left_map[right_id] = max_left;

  return;
}

void EvalLinksV1::calc_max_right(unsigned int left_id) const {
  
  // recalculate max link going right
  std::set<unsigned int> candidates = right(left_id);

  unsigned int max_right = 0xFFFFFFF;
  float max_weight = FLT_MIN;

  for (std::set<unsigned int>::const_iterator iter = candidates.begin();
       iter != candidates.end();
       ++iter) {
    unsigned int candidate = *iter;
    float candidate_weight = get_weight(left_id,candidate);
    
    if (candidate_weight > max_weight) {
      max_right = candidate;
      max_weight = candidate_weight;
    }
  }

  _left_right_map[left_id] = max_right;
}
