
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
const unsigned long long EvalLinksV1::NULLID = ULONG_LONG_MAX;

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
  for (std::map< std::pair<unsigned long long,unsigned long long>, double>::const_iterator citer = _links.begin();
       citer != _links.end();
       ++citer) {
    unsigned long long left_id = citer->first.first;
    unsigned long long right_id = citer->first.second;
    double weight = citer->second;
    
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

void EvalLinksV1::link(unsigned long long left_id, 
		       unsigned long long right_id,
		       double weight) {
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

void EvalLinksV1::unlink(unsigned long long left_id, unsigned long long right_id) {
  if (stale()) refresh();

  if (_links.find(make_pair(left_id,right_id)) != _links.end()) {
    _links.erase(make_pair(left_id,right_id));
    refresh();
  }
  
  return;
}

void EvalLinksV1::unlink_subleading() {
  if (stale()) refresh();

  for (std::map< std::pair<unsigned long long,unsigned long long>, double>::const_iterator citer = _links.begin();
       citer != _links.end();
       ++citer) {
    unsigned long long left_id = citer->first.first;
    unsigned long long right_id = citer->first.second;
    unsigned long long max_right_id = max_right(left_id);

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

bool EvalLinksV1::has_link(unsigned long long left_id, unsigned long long right_id) const {
  if (stale()) refresh();

  if (_links.find(make_pair(left_id,right_id)) != _links.end()) {
    return true;
  }

  return false;
}

double EvalLinksV1::get_weight(unsigned long long left_id, unsigned long long right_id) const {
  if (stale()) refresh();

  std::pair<unsigned long long,unsigned long long> pair = make_pair(left_id,right_id);
  std::map<std::pair<unsigned long long,unsigned long long>, double >::const_iterator citer = _links.find(pair);
  if (citer != _links.end()) {
    return citer->second;
  }

  return NAN;
}

std::set<unsigned long long> EvalLinksV1::left(unsigned long long right_id) const {
  if (stale()) refresh();
  
  std::set<unsigned long long> left_links;
  for (std::multimap<unsigned long long,unsigned long long>::const_iterator citer = _right_left_mmap.lower_bound(right_id);
       citer != _right_left_mmap.upper_bound(right_id);
       ++citer) {
    left_links.insert(citer->second);
  }
  
  return left_links;
}

std::set<unsigned long long> EvalLinksV1::right(unsigned long long left_id) const {
  if (stale()) refresh();
   
  std::set<unsigned long long> right_links;
  for (std::multimap<unsigned long long,unsigned long long>::const_iterator citer = _left_right_mmap.lower_bound(left_id);
       citer != _left_right_mmap.upper_bound(left_id);
       ++citer) {
    right_links.insert(citer->second);
  }
  
  return right_links;
}

unsigned long long EvalLinksV1::max_left(unsigned long long right_id) const {
  if (stale()) refresh();
  if (_right_left_map.find(right_id) == _right_left_map.end()) {
    return NULLID;
  }
  return _right_left_map[right_id];
}

unsigned long long EvalLinksV1::max_right(unsigned long long left_id) const {
  if (stale()) refresh();
  if (_left_right_map.find(left_id) == _left_right_map.end()) {
    return NULLID;
  }
  return _left_right_map[left_id];
}

void EvalLinksV1::refresh() const {

  _left_right_mmap.clear();
  _right_left_mmap.clear();
  _left_right_map.clear();
  _right_left_map.clear();

  // loop over all link pairs and recreate multimaps
  for (std::map<std::pair<unsigned long long,unsigned long long>,double>::const_iterator citer = _links.begin();
       citer != _links.end();
       ++citer) {
    unsigned long long left_id = citer->first.first;
    unsigned long long right_id = citer->first.second;

    _left_right_mmap.insert(make_pair(left_id,right_id));
    _right_left_mmap.insert(make_pair(right_id,left_id));    
  }

  // loop over all multimap keys and recompute max links
  for (std::map<unsigned long long,unsigned long long>::const_iterator citer = _left_right_mmap.begin();
       citer != _left_right_mmap.end();
       ++citer) {
    unsigned long long left_id = citer->first;
    calc_max_right(left_id);
  }

  for (std::map<unsigned long long,unsigned long long>::const_iterator citer = _right_left_mmap.begin();
       citer != _right_left_mmap.end();
       ++citer) {
    unsigned long long right_id = citer->first;
    calc_max_left(right_id);
  }
  
  _stale = false;
}

void EvalLinksV1::calc_max_left(unsigned long long right_id) const {
  
  std::set<unsigned long long> candidates = left(right_id);

  unsigned long long max_left = NULLID;
  double max_weight = DBL_MIN;

  for (std::set<unsigned long long>::const_iterator iter = candidates.begin();
       iter != candidates.end();
       ++iter) {
    unsigned long long candidate = *iter;
    double candidate_weight = get_weight(candidate,right_id);
    
    if (candidate_weight > max_weight) {
      max_left = candidate;
      max_weight = candidate_weight;
    }
  }

  _right_left_map[right_id] = max_left;

  return;
}

void EvalLinksV1::calc_max_right(unsigned long long left_id) const {
  
  // recalculate max link going right
  std::set<unsigned long long> candidates = right(left_id);

  unsigned long long max_right = NULLID;
  double max_weight = DBL_MIN;

  for (std::set<unsigned long long>::const_iterator iter = candidates.begin();
       iter != candidates.end();
       ++iter) {
    unsigned long long candidate = *iter;
    double candidate_weight = get_weight(left_id,candidate);
    
    if (candidate_weight > max_weight) {
      max_right = candidate;
      max_weight = candidate_weight;
    }
  }

  _left_right_map[left_id] = max_right;
}
