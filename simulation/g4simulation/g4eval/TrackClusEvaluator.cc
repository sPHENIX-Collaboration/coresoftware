#include "TrackClusEvaluator.h"
#include "ClusKeyIter.h"
#include "TrkrClusterIsMatcher.h"
#include "g4evalfn.h"

#include <g4tracking/TrkrTruthTrack.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>

#include <phool/phool.h>
#include <numeric>

TrkrClusterContainer* TrackClusEvaluator::get_PHG4_clusters()
{
  if (ismatcher == nullptr)
  {
    return nullptr;
  }
  else
  {
    return ismatcher->m_TruthClusters;
  }
}

TrkrClusterContainer* TrackClusEvaluator::get_SVTX_clusters()
{
  if (ismatcher == nullptr)
  {
    return nullptr;
  }
  else
  {
    return ismatcher->m_RecoClusters;
  }
}

std::array<int, 5> TrackClusEvaluator::cntclus(Vector& keys)
{
  std::array<int, 5> cnt{0, 0, 0, 0, 0};
  for (auto& it : keys)
  {
    cnt[g4evalfn::trklayer_det(it.first)] += 1;
  }
  for (int i = 0; i < 4; ++i)
  {
    cnt[4] += cnt[i];
  }
  return cnt;
}

int TrackClusEvaluator::addClusKeys(SvtxTrack* track)
{
  svtx_keys.clear();
  for (auto ckey : ClusKeyIter(track))
  {
    svtx_keys.push_back({TrkrDefs::getHitSetKeyFromClusKey(ckey), ckey});
  }
  std::sort(svtx_keys.begin(), svtx_keys.end());
  return svtx_keys.size();
}

std::array<int, 5> TrackClusEvaluator::cnt_matchedclus(Vector& keys, std::vector<bool>& matches)
{
  std::array<int, 5> cnt{0, 0, 0, 0, 0};
  if (keys.size() != matches.size())
  {
    std::cout << PHWHERE << " matching and key vector not the same size. "
              << std::endl
              << " run find_matches() first." << std::endl;
    return cnt;
  }
  for (unsigned int i = 0; i < keys.size(); ++i)
  {
    if (matches[i])
    {
      cnt[g4evalfn::trklayer_det(keys[i].first)] += 1;
    }
  }
  for (int i = 0; i < 4; ++i)
  {
    cnt[4] += cnt[i];
  }
  return cnt;
}

int TrackClusEvaluator::addClusKeys(TrkrTruthTrack* track)
{
  phg4_keys.clear();
  for (auto ckey : track->getClusters())
  {
    phg4_keys.push_back({TrkrDefs::getHitSetKeyFromClusKey(ckey), ckey});
  }
  std::sort(phg4_keys.begin(), phg4_keys.end());
  return phg4_keys.size();
}

void TrackClusEvaluator::reset()
{
  phg4_keys.clear();
  phg4_matches.clear();
  svtx_keys.clear();
  svtx_matches.clear();
}

std::array<int, 3> TrackClusEvaluator::find_matches()
{
  if (ismatcher == nullptr)
  {
    std::cout << PHWHERE
              << " Won't compare tracks because of missing TrkrClusterIsMatcher" << std::endl;
    return {0, 0, 0};
  }
  // find the matches between the svtx_keys and phg4_keys
  // also keep track of the sum of the comparison between then

  // ---------------------------------
  // set aliases for notation cleaness
  // use A for PHG4 and B for SVTX
  auto& vA = phg4_keys;
  auto& vB = svtx_keys;

  auto& matchesA = phg4_matches;
  auto& matchesB = svtx_matches;

  match_stat = 0.;  // DEPRECATED

  // matches will say, cluster by cluster, which clusters are matched
  matchesA = std::vector<bool>(vA.size(), false);
  matchesB = std::vector<bool>(vB.size(), false);

  // user iterators to access the vectors
  auto iA0 = vA.begin();
  auto iA1 = vA.end();

  auto iB0 = vB.begin();
  auto iB1 = vB.end();

  auto iA = iA0;
  auto iB = iB0;

  int n_match{0};

  while (iA != iA1 && iB != iB1)
  {
    if (iA->first == iB->first)
    {
      auto hitset = iA->first;

      // must compare ALL sets of iA and iB with this same hitset
      auto sAend = iA + 1;  // search A end
      while (sAend != iA1 && sAend->first == hitset)
      {
        ++sAend;
      }

      auto sBend = iB + 1;  // search B end
      while (sBend != iB1 && sBend->first == hitset)
      {
        ++sBend;
      }

      for (auto A_iter = iA; A_iter != sAend; ++A_iter)
      {
        for (auto B_iter = iB; B_iter != sBend; ++B_iter)
        {
          auto is_match = ismatcher->operator()(A_iter->second, B_iter->second);
          if (is_match)
          {
            matchesA[A_iter - iA0] = true;
            matchesB[B_iter - iB0] = true;
            if (collect_match_statistic)
            {
              match_stat += g4evalfn::calc_match_statistic(ismatcher, A_iter->second, B_iter->second);
            }
            ++n_match;
          }
        }
      }
      iA = sAend;
      iB = sBend;
    }
    else if (iA->first < iB->first)
    {
      ++iA;
    }
    else
    {
      ++iB;
    }
  }
  return {n_match, (int) phg4_keys.size(), (int) svtx_keys.size()};
}

std::array<int, 3> TrackClusEvaluator::find_matches(TrkrTruthTrack* g4_track, SvtxTrack* sv_track)
{
  addClusKeys(sv_track);
  addClusKeys(g4_track);
  return find_matches();
}

int TrackClusEvaluator::phg4_n_matched()
{
  return std::accumulate(phg4_matches.begin(), phg4_matches.end(), 0);
}

int TrackClusEvaluator::svtx_n_matched()
{
  return std::accumulate(svtx_matches.begin(), svtx_matches.end(), 0);
}

std::vector<TrkrClusLoc> TrackClusEvaluator::phg4_clusloc_all()
{
  std::vector<TrkrClusLoc> vec{};
  for (auto& cluspair : phg4_keys)
  {
    vec.push_back(g4evalfn::clusloc_PHG4(ismatcher, cluspair.second));
  }
  return vec;
}

std::vector<TrkrClusLoc> TrackClusEvaluator::phg4_clusloc_unmatched()
{
  std::vector<TrkrClusLoc> vec{};
  auto cnt = phg4_keys.size();
  for (unsigned int i = 0; i < cnt; ++i)
  {
    if (!phg4_matches[i])
    {
      vec.push_back(g4evalfn::clusloc_PHG4(ismatcher, phg4_keys[i].second));
    }
  }
  return vec;
}

std::vector<TrkrClusLoc> TrackClusEvaluator::svtx_clusloc_all()
{
  std::vector<TrkrClusLoc> vec{};
  for (auto& cluspair : svtx_keys)
  {
    vec.push_back(g4evalfn::clusloc_SVTX(ismatcher, cluspair.second));
  }
  return vec;
}

std::vector<TrkrClusLoc> TrackClusEvaluator::svtx_clusloc_unmatched()
{
  std::vector<TrkrClusLoc> vec{};
  auto cnt = svtx_keys.size();
  for (unsigned int i = 0; i < cnt; ++i)
  {
    if (!svtx_matches[i])
    {
      vec.push_back(g4evalfn::clusloc_SVTX(ismatcher, svtx_keys[i].second));
    }
  }
  return vec;
}

std::vector<TrkrClusLoc> TrackClusEvaluator::clusloc_matched()
{
  std::vector<TrkrClusLoc> vec{};
  auto cnt = phg4_keys.size();
  for (unsigned int i = 0; i < cnt; ++i)
  {
    if (phg4_matches[i])
    {
      vec.push_back(g4evalfn::clusloc_PHG4(ismatcher, phg4_keys[i].second));
    }
  }
  return vec;
}
