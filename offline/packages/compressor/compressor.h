/**
 * Compression algorithm for use in the CMSSW and sPHENIX projects.
 * Using 16-bit short integers to store 32-bit floating-point values.
 * Author: fishyu@iii.org.tw
 * May 22, 2021
 */
#include <TTree.h>

#include <fstream>
#include <map>
#include <set>
#include <vector>

//-----------------------------------------------------------------------------
/**
 * approx() compresses data held in t and returns the standard deviation of the differences between the actual and approximated data.
 */
Float_t approx(
  std::vector<UShort_t>* order, 
  std::vector<Float_t>* dict, 
  std::vector<size_t>* cnt, 
  Int_t n_entries, 
  TTree* t, 
  Float_t* gen_, 
  size_t maxNumClusters
);
//-----------------------------------------------------------------------------
Int_t newLoc(std::vector<Int_t>* loc_vec, std::vector<std::set<Int_t>>* loc_set_vec);
void removeDiff(Float_t distance, Float_t min, std::map<Float_t, std::set<Float_t>>* distance_min_set_map);
void addDiff(Float_t distance, Float_t min, std::map<Float_t, std::set<Float_t>>* distance_min_set_map);
//-----------------------------------------------------------------------------
Float_t approx(std::vector<UShort_t>* order, std::vector<Float_t>* dict, std::vector<size_t>* cnt, Int_t n_entries, TTree* t, Float_t* gen_, size_t maxNumClusters)
{
  Float_t maxAbsErrorDoubled = (Float_t) 0;

  std::map<Float_t, std::pair<Float_t, Int_t>> min_max_loc_map;
  std::vector<std::set<Int_t>> loc_set_vec;
  std::vector<Int_t> loc_vec;
  std::map<Float_t, std::set<Float_t>> distance_min_set_map;

  for (Int_t j = 0 ; j < n_entries; j++){
    t->GetEntry(j);
    
    std::map<Float_t, std::pair<Float_t, Int_t>>::iterator mmlm = min_max_loc_map.find(*gen_);

    if (mmlm != min_max_loc_map.end())
      loc_set_vec[mmlm->second.second].insert(j);
    else {
      Int_t loc = newLoc(&loc_vec, &loc_set_vec);

      loc_set_vec[loc].insert(j);

      min_max_loc_map[*gen_] = std::pair<Float_t, Int_t>(*gen_, loc);

      mmlm = min_max_loc_map.find(*gen_);
      if (mmlm != min_max_loc_map.begin() && *gen_ <= prev(mmlm)->second.first) {
        loc_set_vec[prev(mmlm)->second.second].insert(j);
        loc_set_vec[mmlm->second.second].clear();
        loc_vec.push_back(mmlm->second.second);

        min_max_loc_map.erase(mmlm);

      } else if (min_max_loc_map.size() >= 2) {
        if (mmlm != min_max_loc_map.begin() && mmlm != prev(min_max_loc_map.end())) {

          removeDiff(next(mmlm)->second.first - prev(mmlm)->first, prev(mmlm)->first, &distance_min_set_map);
        }

        if (mmlm != min_max_loc_map.begin())
          addDiff(mmlm->second.first - prev(mmlm)->first, prev(mmlm)->first, &distance_min_set_map);

        if (mmlm != prev(min_max_loc_map.end()))
          addDiff(next(mmlm)->second.first - mmlm->first, mmlm->first, &distance_min_set_map);
      }
    }

    if (min_max_loc_map.size() <= maxNumClusters)
      continue;
   
    std::map<Float_t, std::set<Float_t>>::iterator dmsm = distance_min_set_map.begin();
    Float_t min = *(dmsm->second.begin());
    
    dmsm->second.erase(min);
    if (dmsm->second.empty())
      distance_min_set_map.erase(dmsm);
    
    mmlm = min_max_loc_map.find(min);
    if (mmlm != min_max_loc_map.begin())
      removeDiff(mmlm->second.first - prev(mmlm)->first, prev(mmlm)->first, &distance_min_set_map);

    if (next(mmlm) != prev(min_max_loc_map.end()))
      removeDiff(next(next(mmlm))->second.first - next(mmlm)->first, next(mmlm)->first, &distance_min_set_map);

    std::set<Int_t>* s = &(loc_set_vec[next(mmlm)->second.second]);
    loc_set_vec[mmlm->second.second].insert(s->begin(), s->end());
    mmlm->second.first = next(mmlm)->second.first;
    min_max_loc_map.erase(next(mmlm));
    mmlm = min_max_loc_map.find(min);
    maxAbsErrorDoubled = std::max(maxAbsErrorDoubled, mmlm->second.first - mmlm->first);
    if (mmlm != min_max_loc_map.begin())
      addDiff(mmlm->second.first - prev(mmlm)->first, prev(mmlm)->first, &distance_min_set_map);

    if (mmlm != prev(min_max_loc_map.end()))
      addDiff(next(mmlm)->second.first - mmlm->first, mmlm->first, &distance_min_set_map);

  }
  
  Double_t squaredSum = 0;
  Double_t sum = 0;
  
  order->resize(n_entries);
  for (const auto &mmlm : min_max_loc_map) {
	Double_t estimate = (Double_t) (mmlm.first + mmlm.second.first) / (Double_t) 2;
	  
    for (const auto &index : loc_set_vec[mmlm.second.second]) {
      (*order)[index] = dict->size();
	  
	  t->GetEntry(index);
	  Double_t delta = std::fabs(*gen_ - estimate);
	  squaredSum += (delta * delta);
	  sum += delta;
	}

    dict->push_back(estimate);
    cnt->push_back(loc_set_vec[mmlm.second.second].size());
  }
  
  Double_t avg = sum / (Double_t) n_entries;
  return sqrt((squaredSum / (Double_t) n_entries) - avg * avg);
}

Int_t newLoc(std::vector<Int_t>* loc_vec, std::vector<std::set<Int_t>>* loc_set_vec)
{
  if (!loc_vec->empty()) {
    Int_t loc = loc_vec->back();
    loc_vec->pop_back();
    return loc;
  }
 
  Int_t loc = loc_set_vec->size();
  loc_set_vec->push_back({});
  return loc;
}


void removeDiff(Float_t distance, Float_t min, std::map<Float_t, std::set<Float_t>>* distance_min_set_map) 
{
  std::map<Float_t, std::set<Float_t>>::iterator dmsm = distance_min_set_map->find(distance);
  dmsm->second.erase(min);

  if (dmsm->second.empty())
    distance_min_set_map->erase(dmsm);
}

void addDiff(Float_t distance, Float_t min, std::map<Float_t, std::set<Float_t>>* distance_min_set_map) 
{
  std::map<Float_t, std::set<Float_t>>::iterator dmsm = distance_min_set_map->find(distance);
  if (dmsm == distance_min_set_map->end()) {
    (*distance_min_set_map)[distance] = {min};
  } else {
    dmsm->second.insert(min);
  }
}
