/*!
 *  \file PHCASeeding.cc
 *  \brief Track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail
 *  \author Michael Peters & Christof Roland
 */

#include "PHCASeeding.h"
#include "ALICEKF.h"
#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHTimer.h>  // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

// tpc distortion correction
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <phfield/PHFieldConfigv2.h>

#include <ffamodules/CDBInterface.h>

// trackbase_historic includes
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed_v2.h>

// ROOT includes for debugging
#include <TFile.h>
#include <TNtuple.h>

// BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/policies/compare.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <algorithm>
#include <algorithm>  // for find
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <numeric>
#include <unordered_set>
#include <utility>  // for pair, make_pair
#include <vector>


//#define _DEBUG_

#if defined(_DEBUG_)
#define LogDebug(exp) \
  if (Verbosity() > 0) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void) 0
#endif

#define LogError(exp) \
  if (Verbosity() > 0) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) \
  if (Verbosity() > 0) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

// _CLUSTER_LOG_TUPOUT_ defined statement in the header file
#if defined(_CLUSTER_LOG_TUPOUT_)
#define _FILL_TUPLE(tupname, num, key, pos) \
  tupname->Fill(_nevent,TrkrDefs::getLayer(key), num, pos.x(), pos.y(), pos.z())
#define _FILL_TUPLE_WITH_SEED(tupname, seed,pos) \
  for (unsigned int i=0; i<seed.size();++i) _FILL_TUPLE(tupname, i, seed[i], pos.at(seed[i]))
#define _PROGRESS_TUPOUT_COUNT() _tupout_count += 1
#else 
#define _FILL_TUPLE(tupname, num, key,pos) (void) 0
#define _FILL_TUPLE_WITH_SEED(tupname, seed, pos) (void) 0
#define _PROGRESS_TUPOUT_COUNT() (void) 0
#endif

#if defined(_PHCASEEDING_TIMER_OUT_)
#define _PHCASEEDING_PRINT_TIME(timer,statement) timer->stop(); \
  std::cout << " _PHCASEEDING_PRINT_TIME: Time to " << statement << ": " << timer->elapsed()/1000 << " s" << std::endl;
#else
#define _PHCASEEDING_PRINT_TIME(timer, statement) (void) 0
#endif

//#define _DEBUG_

// end


// apparently there is no builtin STL hash function for a std::array
// so to use std::unordered_set (essentially a hash table), we have to make our own hasher

namespace std
{
  template <typename T, size_t N>
  struct hash<std::array<T, N>>
  {
    using argument_type = std::array<T, N>;
    using result_type = size_t;

    result_type operator()(const argument_type& a) const
    {
      hash<T> hasher;
      result_type h = 0;
      for (result_type i = 0; i < N; ++i)
      {
        h = h * 31 + hasher(a[i]);
      }
      return h;
    }
  };
  template <typename A, typename B>
  struct hash<pair<A, B>>
  {
    using argument_type = pair<A, B>;
    using result_type = size_t;

    result_type operator()(const argument_type& a) const
    {
      hash<A> hashA;
      hash<B> hashB;
      return (hashA(a.first) * 31 + hashB(a.second));
    }
  };
}  // namespace std

// anonymous namespace for local functions
namespace
{
  // square
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  /// phi angle of Acts::Vector3
  inline double get_phi(const Acts::Vector3& position)
  {
    double phi = std::atan2(position.y(), position.x());
    if (phi < 0)
    {
      phi += 2. * M_PI;
    }
    return phi;
  }

  /// pseudo rapidity of Acts::Vector3
  inline double get_eta(const Acts::Vector3& position)
  {
    const double norm = std::sqrt(square(position.x()) + square(position.y()) + square(position.z()));
    return std::log((norm + position.z()) / (norm - position.z())) / 2;
  }

  inline double breaking_angle(double x1, double y1, double z1, double x2, double y2, double z2)
  {
    double l1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double l2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    double sx = (x1 / l1 + x2 / l2);
    double sy = (y1 / l1 + y2 / l2);
    double sz = (z1 / l1 + z2 / l2);
    double dx = (x1 / l1 - x2 / l2);
    double dy = (y1 / l1 - y2 / l2);
    double dz = (z1 / l1 - z2 / l2);
    return 2 * atan2(sqrt(dx * dx + dy * dy + dz * dz), sqrt(sx * sx + sy * sy + sz * sz));
  }

}  // namespace

// using namespace ROOT::Minuit2;
namespace bgi = boost::geometry::index;

PHCASeeding::PHCASeeding(
    const std::string& name,
    unsigned int start_layer,
    unsigned int end_layer,
    unsigned int min_nhits_per_cluster,
    unsigned int min_clusters_per_track,
    unsigned int nlayers_maps,
    unsigned int nlayers_intt,
    unsigned int nlayers_tpc,
    float neighbor_phi_width,
    float neighbor_eta_width,
    float maxSinPhi,
    float cosTheta_limit)
  : PHTrackSeeding(name)
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _start_layer(start_layer)
  , _end_layer(end_layer)
  , _min_nhits_per_cluster(min_nhits_per_cluster)
  , _min_clusters_per_track(min_clusters_per_track)
  , _neighbor_phi_width(neighbor_phi_width)
  , _neighbor_eta_width(neighbor_eta_width)
  , _max_sin_phi(maxSinPhi)
  , _cosTheta_limit(cosTheta_limit)
{
}

int PHCASeeding::InitializeGeometry(PHCompositeNode* topNode)
{
  tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cout << PHWHERE << "No acts tracking geometry, can't proceed" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 PHCASeeding::getGlobalPosition(TrkrDefs::cluskey key, TrkrCluster* cluster) const
{
  // get global position from Acts transform
  auto globalpos = tGeometry->getGlobalPosition(key, cluster);

  // check if TPC distortion correction are in place and apply
  if (m_dcc && !(_pp_mode))
  {
    globalpos = m_distortionCorrection.get_corrected_position(globalpos, m_dcc);
  }

  return globalpos;
}

void PHCASeeding::QueryTree(const bgi::rtree<PHCASeeding::pointKey, bgi::quadratic<16>>& rtree, double phimin, double etamin, double phimax, double etamax, std::vector<pointKey>& returned_values) const
{
  bool query_both_ends = false;
  if (phimin < 0)
  {
    query_both_ends = true;
    phimin += 2 * M_PI;
  }
  if (phimax > 2 * M_PI)
  {
    query_both_ends = true;
    phimax -= 2 * M_PI;
  }
  if (query_both_ends)
  {
    rtree.query(bgi::intersects(box(point(phimin, etamin), point(2 * M_PI, etamax))), std::back_inserter(returned_values));
    rtree.query(bgi::intersects(box(point(0., etamin), point(phimax, etamax))), std::back_inserter(returned_values));
  }
  else
  {
    rtree.query(bgi::intersects(box(point(phimin, etamin), point(phimax, etamax))), std::back_inserter(returned_values));
  }
}

std::pair<PHCASeeding::PositionMap, PHCASeeding::keyListPerLayer> PHCASeeding::FillGlobalPositions() {
  keyListPerLayer ckeys;
  PositionMap cachedPositions;
  cachedPositions.reserve(_cluster_map->size());  //avoid resizing mid-execution

  for (const auto& hitsetkey : _cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto range = _cluster_map->getClusters(hitsetkey);
    for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
    {
      TrkrDefs::cluskey ckey = clusIter->first;
      TrkrCluster* cluster = clusIter->second;
      unsigned int layer = TrkrDefs::getLayer(ckey);
      if (layer < _start_layer || layer >= _end_layer)
      {
        if (Verbosity() > 0)
        {
          std::cout << "layer: " << layer << std::endl;
        }
        continue;
      }
      if (_iteration_map != nullptr && _n_iteration > 0)
      {
        if (_iteration_map->getIteration(ckey) > 0)
        {
          continue;  // skip hits used in a previous iteration
        }
      }

      // get global position, convert to Acts::Vector3 and store in map
      const Acts::Vector3 globalpos_d = getGlobalPosition(ckey, cluster);
      const Acts::Vector3 globalpos = {globalpos_d.x(), globalpos_d.y(), globalpos_d.z()};
      cachedPositions.insert(std::make_pair(ckey, globalpos));
      ckeys[layer-_FIRST_LAYER_TPC].push_back(ckey);
      _FILL_TUPLE(_tupclus_all, 0, ckey, globalpos);
    }
  }
  return std::make_pair(cachedPositions, ckeys);
}

std::vector<PHCASeeding::coordKey> PHCASeeding::FillTree(bgi::rtree<PHCASeeding::pointKey,bgi::quadratic<16>>& _rtree, const PHCASeeding::keyList& ckeys, const PHCASeeding::PositionMap& globalPositions, const int layer) 
{
  // Fill _rtree with the clusters in ckeys; remove duplicates, and return a vector of the coordKeys
  // Note that layer is only used for a cout statement
  int n_dupli = 0;
  std::vector<coordKey> coords;
  _rtree.clear();
  /* _rtree.reserve(ckeys.size()); */
  for (const auto& ckey : ckeys) {
    const auto& globalpos_d = globalPositions.at(ckey);
    const double clus_phi = get_phi(globalpos_d);
    const double clus_eta = get_eta(globalpos_d);
    if (Verbosity() > 0)
    {
      /* int layer = TrkrDefs::getLayer(ckey); */
      std::cout << "Found cluster " << ckey << " in layer " << layer << std::endl;
    }
    std::vector<pointKey> testduplicate;
    QueryTree(_rtree, clus_phi - 0.00001, clus_eta - 0.00001, clus_phi + 0.00001, clus_eta + 0.00001, testduplicate);
    if (!testduplicate.empty())
    {
      ++n_dupli;
      continue;
    }
    coords.push_back({{static_cast<float>(clus_phi), static_cast<float>(clus_eta)},ckey});
    t_fill->restart();
    _rtree.insert(std::make_pair(point(clus_phi, globalpos_d.z()), ckey));
    t_fill->stop();
  }
  if (Verbosity() > 1)
  {
      std::cout << "nhits in layer("<<layer<<"): " << coords.size() << std::endl;
  }
  if (Verbosity() > 0)
  {
    std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
  }
  if (Verbosity() > 0)
  {
    std::cout << "number of duplicates : " << n_dupli << std::endl;
  }
  return coords;
}

int PHCASeeding::Process(PHCompositeNode* /*topNode*/)
{
  _PROGRESS_TUPOUT_COUNT();
  if (Verbosity() > 1)
  {
    std::cout << " Process...  " << std::endl;
  }
  //  TFile fpara("CA_para.root", "RECREATE");
  if (_n_iteration > 0)
  {
    if (!_iteration_map)
    {
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  t_seed->restart();
  t_makebilinks->restart();

  PositionMap globalPositions;
  keyListPerLayer ckeys;
  std::tie(globalPositions, ckeys) = FillGlobalPositions();

  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "Initial RTree fill time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  t_seed->restart();
  int numberofseeds = 0;
  numberofseeds += FindSeedsWithMerger(globalPositions,ckeys);
  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "number of seeds " << numberofseeds << std::endl;
  }
  if (Verbosity() > 0)
  {
    std::cout << "Kalman filtering time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  //  fpara.cd();
  //  fpara.Close();
  //  if(Verbosity()>0) std::cout << "fpara OK\n";
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::FindSeedsWithMerger(const PHCASeeding::PositionMap& globalPositions, const PHCASeeding::keyListPerLayer& ckeys)
{
  t_seed->restart();

  keyLinks trackSeedPairs;
  keyLinkPerLayer bodyLinks;
  std::tie(trackSeedPairs, bodyLinks) = CreateBiLinks(globalPositions, ckeys);
  _PHCASEEDING_PRINT_TIME(t_makebilinks, "init and make bilinks");

  t_makeseeds->restart();
  keyLists trackSeedKeyLists = FollowBiLinks(trackSeedPairs, bodyLinks, globalPositions);
  _PHCASEEDING_PRINT_TIME(t_makeseeds, "make seeds");
  if (Verbosity() > 0 )
  {
    t_makeseeds->stop();
    std::cout << "Time to make seeds: " << t_makeseeds->elapsed() / 1000 << " s" << std::endl;
  }
  std::vector<TrackSeed_v2> seeds = RemoveBadClusters(trackSeedKeyLists, globalPositions);

  publishSeeds(seeds);
  return seeds.size();

}

std::pair<PHCASeeding::keyLinks, PHCASeeding::keyLinkPerLayer>  PHCASeeding::CreateBiLinks(const PHCASeeding::PositionMap& globalPositions, const PHCASeeding::keyListPerLayer& ckeys)
{
  keyLinks startLinks; // bilinks at start of chains
  keyLinkPerLayer bodyLinks; //  bilinks to build chains
                              //
  double cluster_find_time = 0;
  double rtree_query_time = 0;
  double transform_time = 0;
  double compute_best_angle_time = 0;
  double set_insert_time = 0;


  // there are three coord_array (only the current layer is used at a time,
  // but it is filled the same time as the _rtrees, which are used two at
  // a time -- the prior padplane row and the next padplain row
  std::array<std::vector<coordKey>,3> coord_arr;
  std::array<std::unordered_set<keyLink>,2> previous_downlinks_arr;
  std::array<std::unordered_set<TrkrDefs::cluskey>,2> bottom_of_bilink_arr;

  // iterate from outer to inner layers
  const int inner_index = _start_layer - _FIRST_LAYER_TPC+1;
  const int outer_index = _end_layer   - _FIRST_LAYER_TPC-2;

  // fill the current and prior row coord and ttrees for the first iteration
  int _index_above   = (outer_index+1)%3;
  int _index_current = (outer_index)%3;
  coord_arr[_index_above]   = FillTree(_rtrees[_index_above],   ckeys[outer_index+1], globalPositions, outer_index+1);
  coord_arr[_index_current] = FillTree(_rtrees[_index_current], ckeys[outer_index],   globalPositions, outer_index);

  for (int layer_index=outer_index; layer_index>=inner_index;--layer_index) {
    // these lines of code will rotates through all three _rtree's in the array,
    // where the old upper becomes the new middle, the old middle the new lower,
    // and the old lower drops out and that _rtree is filled with the new upper
    int index_above   = (layer_index+1)%3;
    int index_current = (layer_index  )%3;
    int index_below   = (layer_index-1)%3;

    coord_arr[index_below] = FillTree(_rtrees[index_below], ckeys[layer_index-1], globalPositions, layer_index-1);

    auto& _rtree_above = _rtrees[index_above];
    const std::vector<coordKey>& coord        = coord_arr[index_current];
    auto& _rtree_below = _rtrees[index_below];

    auto& curr_downlinks = previous_downlinks_arr[layer_index%2];
    auto& last_downlinks = previous_downlinks_arr[(layer_index+1)%2];

    auto& curr_bottom_of_bilink = bottom_of_bilink_arr[layer_index%2];
    auto& last_bottom_of_bilink = bottom_of_bilink_arr[(layer_index+1)%2];

    curr_downlinks.clear();
    curr_bottom_of_bilink.clear();

    // For all the clusters in coord, find nearest neighbors in the 
    // above and below layers and make links
    // Any link to an above node which matches the same clusters
    // on the previous iteration (to a "below node") becomes a "bilink"
    // Each bilink will either add to an existing chain or start a new one
    /* std::vector<keyLink> belowLinks; */
    std::vector<keyLink> aboveLinks;
    for (const auto& StartCluster : coord)
    {
      double StartPhi = StartCluster.first[0];
      const auto& globalpos = globalPositions.at(StartCluster.second);
      double StartX = globalpos(0);
      double StartY = globalpos(1);
      double StartZ = globalpos(2);
      t_seed->stop();
      cluster_find_time += t_seed->elapsed();
      t_seed->restart();
      LogDebug(" starting cluster:" << std::endl);
      LogDebug(" z: " << StartZ << std::endl);
      LogDebug(" phi: " << StartPhi << std::endl);

      std::vector<pointKey> ClustersAbove;
      std::vector<pointKey> ClustersBelow;

      QueryTree(_rtree_below,
          StartPhi - _neighbor_phi_width,
          StartZ - _neighbor_eta_width,
          StartPhi + _neighbor_phi_width,
          StartZ + _neighbor_eta_width,
          ClustersBelow);

      QueryTree(_rtree_above,
          StartPhi - _neighbor_phi_width,
          StartZ - _neighbor_eta_width,
          StartPhi + _neighbor_phi_width,
          StartZ + _neighbor_eta_width,
          ClustersAbove);

      t_seed->stop();
      rtree_query_time += t_seed->elapsed();
      t_seed->restart();
      LogDebug(" entries in below layer: " << ClustersBelow.size() << std::endl);
      LogDebug(" entries in above layer: " << ClustersAbove.size() << std::endl);
      std::vector<std::array<double, 3>> delta_below;
      std::vector<std::array<double, 3>> delta_above;
      delta_below.clear();
      delta_above.clear();
      delta_below.resize(ClustersBelow.size());
      delta_above.resize(ClustersAbove.size());
      // calculate (delta_eta, delta_phi) vector for each neighboring cluster

      std::transform(ClustersBelow.begin(), ClustersBelow.end(), delta_below.begin(),
          [&](pointKey BelowCandidate)
          {
          const auto& belowpos = globalPositions.at(BelowCandidate.second);
          return std::array<double,3>{belowpos(0)-StartX,
          belowpos(1)-StartY,
          belowpos(2)-StartZ}; });

      std::transform(ClustersAbove.begin(), ClustersAbove.end(), delta_above.begin(),
          [&](pointKey AboveCandidate)
          {
          const auto& abovepos = globalPositions.at(AboveCandidate.second);
          return std::array<double,3>{abovepos(0)-StartX,
          abovepos(1)-StartY,
          abovepos(2)-StartZ}; });
      t_seed->stop();
      transform_time += t_seed->elapsed();
      t_seed->restart();

      // find the three clusters closest to a straight line
      // (by maximizing the cos of the angle between the (delta_eta,delta_phi) vectors)
      double maxCosPlaneAngle = -0.95;
      // double minSumLengths = 1e9;
      keyList bestBelowClusters;
      keyList bestAboveClusters;
      for (size_t iAbove = 0; iAbove < delta_above.size(); ++iAbove)
      {
        for (size_t iBelow = 0; iBelow < delta_below.size(); ++iBelow)
        {
          //        double dotProduct = delta_below[iBelow][0]*delta_above[iAbove][0]+delta_below[iBelow][1]*delta_above[iAbove][1]+delta_below[iBelow][2]*delta_above[iAbove][2];
          //        double belowLength = sqrt(delta_below[iBelow][0]*delta_below[iBelow][0]+delta_below[iBelow][1]*delta_below[iBelow][1]+delta_below[iBelow][2]*delta_below[iBelow][2]);
          //        double aboveLength = sqrt(delta_above[iAbove][0]*delta_above[iAbove][0]+delta_above[iAbove][1]*delta_above[iAbove][1]+delta_above[iAbove][2]*delta_above[iAbove][2]);
          double angle = breaking_angle(
              delta_below[iBelow][0],
              delta_below[iBelow][1],
              delta_below[iBelow][2],
              delta_above[iAbove][0],
              delta_above[iAbove][1],
              delta_above[iAbove][2]);
          if (cos(angle) < maxCosPlaneAngle)
          {
            // maxCosPlaneAngle = cos(angle);
            // minSumLengths = belowLength+aboveLength;
            bestBelowClusters.push_back(ClustersBelow[iBelow].second);
            bestAboveClusters.push_back(ClustersAbove[iAbove].second);
          }
        }
      }
      // NOTE:
      // There was some old commented-out code here for allowing layers to be skipped. This 
      // may be useful in the future. This chunk of code has been moved towards the 
      // end fo the file under the title: "---OLD CODE 0: SKIP_LAYERS---"
      t_seed->stop();
      compute_best_angle_time += t_seed->elapsed();
      t_seed->restart();

      for (auto cluster : bestBelowClusters)
      {
        curr_downlinks.insert({StartCluster.second, cluster});
        _FILL_TUPLE(_tupclus_links, 0, StartCluster.second, globalPositions.at(cluster));
        _FILL_TUPLE(_tupclus_links, -1, cluster, globalPositions.at(cluster));
      }

      for (auto cluster : bestAboveClusters)
      {
        _FILL_TUPLE(_tupclus_links, 1, cluster, globalPositions.at(cluster));
        keyLink uplink = std::make_pair(cluster, StartCluster.second);

        if (last_downlinks.find(uplink) != last_downlinks.end())
        {
          // this is a bilink
          const auto& key_top = uplink.first;
          const auto& key_bot = uplink.second;
          curr_bottom_of_bilink.insert(key_bot);
          _FILL_TUPLE(_tupclus_bilinks, 0, key_top, globalPositions.at(cluster));
          _FILL_TUPLE(_tupclus_bilinks, 1, key_bot, globalPositions.at(cluster));
          
          if (last_bottom_of_bilink.find(key_top)==last_bottom_of_bilink.end()) {
            startLinks.push_back(std::make_pair(key_top,key_bot));
          } else {
            bodyLinks[layer_index+1].push_back(std::make_pair(key_top,key_bot));
          }
        }
      } // end loop over all up-links
    } // end loop over start clusters
    t_seed->stop();
    set_insert_time += t_seed->elapsed();
    t_seed->restart();
    LogDebug(" max collinearity: " << maxCosPlaneAngle << std::endl);
  } // end loop over layers (to make links)

  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "triplet forming time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
    std::cout << "starting cluster setup: " << cluster_find_time / 1000 << " s" << std::endl;
    std::cout << "RTree query: " << rtree_query_time / 1000 << " s" << std::endl;
    std::cout << "Transform: " << transform_time / 1000 << " s" << std::endl;
    std::cout << "Compute best triplet: " << compute_best_angle_time / 1000 << " s" << std::endl;
    std::cout << "Set insert: " << set_insert_time / 1000 << " s" << std::endl;
  }
  t_seed->restart();

  // sort the body links per layer so that links can be binary-searched per layer
  /* for (auto& layer : bodyLinks) { std::sort(layer.begin(), layer.end()); } */
  return std::make_pair(startLinks, bodyLinks);
}

double PHCASeeding::getMengerCurvature(TrkrDefs::cluskey a, TrkrDefs::cluskey b, TrkrDefs::cluskey c, const PHCASeeding::PositionMap& globalPositions) const
{
  // Menger curvature = 1/R for circumcircle of triangle formed by most recent three clusters
  // We use here 1/R = 2*sin(breaking angle)/(hypotenuse of triangle)
  auto& a_pos = globalPositions.at(a);
  auto& b_pos = globalPositions.at(b);
  auto& c_pos = globalPositions.at(c);
  double hypot_length = sqrt(square<double>(c_pos.x() - a_pos.x()) + square<double>(c_pos.y() - a_pos.y()) + square<double>(c_pos.z() - a_pos.z()));
  double break_angle = breaking_angle(
      a_pos.x() - b_pos.x(),
      a_pos.y() - b_pos.y(),
      a_pos.z() - b_pos.z(),
      c_pos.x() - b_pos.x(),
      c_pos.y() - b_pos.y(),
      c_pos.z() - b_pos.z());
  return 2 * sin(break_angle) / hypot_length;
}

PHCASeeding::keyLists PHCASeeding::FollowBiLinks( const PHCASeeding::keyLinks& trackSeedPairs, const PHCASeeding::keyLinkPerLayer& bilinks, const PHCASeeding::PositionMap& globalPositions) const
{

  // form all possible starting 3-cluster tracks (we need that to calculate curvature)
  std::vector<keyList> trackSeedKeyLists;
  for (auto& startLink : trackSeedPairs)
  {
    TrkrDefs::cluskey trackHead = startLink.second;
    unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - _FIRST_LAYER_TPC;
    // the following call with get iterators to all bilinks which match the head
    for (const auto& matchlink : bilinks[trackHead_layer]) {
    /* auto matched_links = std::equal_range(bilinks[trackHead_layer].begin(), bilinks[trackHead_layer].end(), trackHead, CompKeyToBilink()); */
    /* for (auto matchlink = matched_links.first; matchlink != matched_links.second; ++matchlink) */
    /* { */
      if (matchlink.first != trackHead) continue;
      keyList trackSeedTriplet;
      trackSeedTriplet.push_back(startLink.first);
      trackSeedTriplet.push_back(startLink.second);
      trackSeedTriplet.push_back(matchlink.second);
      trackSeedKeyLists.push_back(trackSeedTriplet);

      _FILL_TUPLE(_tupclus_seeds, 0, startLink.first, globalPositions.at(startLink.first));
      _FILL_TUPLE(_tupclus_seeds, 1, startLink.second, globalPositions.at(startLink.second));
      _FILL_TUPLE(_tupclus_seeds, 2, matchlinkg.second, globalPositions.at(matchlinkg.second));
    }
  }

  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "starting cluster finding time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  t_seed->restart();
  // assemble track cluster chains from starting cluster keys (ordered from outside in)

  // std::cout << "STARTING SEED ASSEMBLY" << std::endl;

  //  std::vector<keyList> trackSeedKeyLists;
  std::vector<keyList> tempSeedKeyLists = trackSeedKeyLists;
  trackSeedKeyLists.clear();

  while (tempSeedKeyLists.size() > 0)
  {
    if (Verbosity() > 0)
    {
      std::cout << "temp size: " << tempSeedKeyLists.size() << std::endl;
    }
    if (Verbosity() > 0)
    {
      std::cout << "final size: " << trackSeedKeyLists.size() << std::endl;
    }
    std::vector<keyList> newtempSeedKeyLists;
    for (auto& seed : tempSeedKeyLists)
    {
      TrkrDefs::cluskey trackHead = seed.back();
      unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - (_nlayers_intt + _nlayers_maps);
      // bool no_next_link = true;
      /* for (auto testlink = matched_links.first; testlink != matched_links.second; ++testlink) */
    for (const auto& link : bilinks[trackHead_layer]) {
      if (link.first != trackHead) continue;
      // It appears that it is just faster to traverse the lists, then use a binary-sorted search
      // In any case, if we use this cord in the future, be sure to sort the bilinks before using
    /* auto matched_links = std::equal_range(bilinks[trackHead_layer].begin(), bilinks[trackHead_layer].end(), trackHead, CompKeyToBilink()); */
    /* for (auto link = matched_links.first; link != matched_links.second; ++link) */
      /* { */
        auto& head_pos = globalPositions.at(trackHead);
        auto& prev_pos = globalPositions.at(seed.rbegin()[1]);
        float x1 = head_pos.x();
        float y1 = head_pos.y();
        float z1 = head_pos.z();
        float x2 = prev_pos.x();
        float y2 = prev_pos.y();
        float z2 = prev_pos.z();
        float dr_12 = sqrt(x1 * x1 + y1 * y1) - sqrt(x2 * x2 + y2 * y2);
        TrkrDefs::cluskey testCluster = link.second;
        auto& test_pos = globalPositions.at(testCluster);
        float xt = test_pos.x();
        float yt = test_pos.y();
        float zt = test_pos.z();
        float new_dr = sqrt(xt * xt + yt * yt) - sqrt(x1 * x1 + y1 * y1);
        if (fabs((z1 - z2) / dr_12 - (zt - z1) / new_dr) > _clusadd_delta_dzdr_window)
        { continue; }
        auto& third_pos = globalPositions.at(seed.rbegin()[2]);
        float x3 = third_pos.x();
        float y3 = third_pos.y();
        float dr_23 = sqrt(x2 * x2 + y2 * y2) - sqrt(x3 * x3 + y3 * y3);
        float phi1 = atan2(y1, x1);
        float phi2 = atan2(y2, x2);
        float phi3 = atan2(y3, x3);
        float dphi12 = std::fmod(phi1 - phi2, M_PI);
        float dphi23 = std::fmod(phi2 - phi3, M_PI);
        float d2phidr2 = dphi12 / (dr_12 * dr_12) - dphi23 / (dr_23 * dr_23);
        float new_dphi = std::fmod(atan2(yt, xt) - atan2(y1, x1), M_PI);
        float new_d2phidr2 = new_dphi / (new_dr * new_dr) - dphi12 / (dr_12 * dr_12);
        if (seed.size() < 6 && fabs(d2phidr2 - new_d2phidr2) < _clusadd_delta_dphidr2_window)
        {
          //    no_next_link = false;
          keyList newseed = seed;
          newseed.push_back(link.second);
          newtempSeedKeyLists.push_back(newseed);
        }
      }
      if (seed.size() > 5)
      {
        trackSeedKeyLists.push_back(seed);
        _FILL_TUPLE_WITH_SEED(_tupclus_grown_seeds, seed, globalPositions);
      }
    }
    if (Verbosity() > 0)
    {
      std::cout << "new temp size: " << newtempSeedKeyLists.size() << std::endl;
    }
    tempSeedKeyLists = newtempSeedKeyLists;
  }
  // old code block move to end of code under the title: "---OLD CODE 1: SKIP_LAYERS---"
  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "keychain assembly time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  t_seed->restart();
  LogDebug(" track key chains assembled: " << trackSeedKeyLists.size() << std::endl);
  LogDebug(" track key chain lengths: " << std::endl);
  for (auto trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    LogDebug(" " << trackKeyChain->size() << std::endl);
  }
  int jumpcount = 0;
  LogDebug(" track key associations:" << std::endl);
  for (auto& trackSeedKeyList : trackSeedKeyLists)
  {
    LogDebug(" seed " << i << ":" << std::endl);

    double lasteta = -100;
    double lastphi = -100;
    for (unsigned long& j : trackSeedKeyList)
    {
      const auto& globalpos = globalPositions.at(j); 
      const double clus_phi = get_phi(globalpos);
      const double clus_eta = get_eta(globalpos);
      const double etajump = clus_eta - lasteta;
      const double phijump = clus_phi - lastphi;
#if defined(_DEBUG_)
      unsigned int lay = TrkrDefs::getLayer(trackSeedKeyLists[i][j].second);
#endif
      if ((fabs(etajump) > 0.1 && lasteta != -100) || (fabs(phijump) > 1 && lastphi != -100))
      {
        LogDebug(" Eta or Phi jump too large! " << std::endl);
        ++jumpcount;
      }
      LogDebug(" (eta,phi,layer) = (" << clus_eta << "," << clus_phi << "," << lay << ") "
                                      << " (x,y,z) = (" << globalpos(0) << "," << globalpos(1) << "," << globalpos(2) << ")" << std::endl);

      if (Verbosity() > 0)
      {
        unsigned int lay = TrkrDefs::getLayer(j);
        std::cout << "  eta, phi, layer = (" << clus_eta << "," << clus_phi << "," << lay << ") "
                  << " (x,y,z) = (" << globalpos(0) << "," << globalpos(1) << "," << globalpos(2) << ")" << std::endl;
      }
      lasteta = clus_eta;
      lastphi = clus_phi;
    }
  }
  LogDebug(" Total large jumps: " << jumpcount << std::endl);
  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "eta-phi sanity check time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  t_seed->restart();
  return trackSeedKeyLists;
}

std::vector<TrackSeed_v2> PHCASeeding::RemoveBadClusters(const std::vector<PHCASeeding::keyList>& chains, const PHCASeeding::PositionMap& globalPositions) const
{
  if (Verbosity() > 0)
  {
    std::cout << "removing bad clusters" << std::endl;
  }
  std::vector<TrackSeed_v2> clean_chains;

  for (const auto& chain : chains)
  {
    if (chain.size() < 3)
    {
      continue;
    }
    if (Verbosity() > 0)
    {
      std::cout << "chain size: " << chain.size() << std::endl;
    }

    TrackFitUtils::position_vector_t xy_pts;
    for (const auto& cluskey : chain)
    {
      const auto& global = globalPositions.at(cluskey);
      xy_pts.emplace_back(global.x(), global.y());
    }

    // fit a circle through x,y coordinates
    const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin(xy_pts);

    // skip chain entirely if fit fails
    if (std::isnan(R))
    {
      continue;
    }

    // calculate residuals
    const std::vector<double> xy_resid = TrackFitUtils::getCircleClusterResiduals(xy_pts, R, X0, Y0);

    // assign clusters to seed
    TrackSeed_v2 trackseed;
    for (unsigned long i : chain)
    {
      // if(xy_resid[i]>_xy_outlier_threshold) continue;
      trackseed.insert_cluster_key(i);
    }

    clean_chains.push_back(trackseed);
    if (Verbosity() > 0)
    {
      std::cout << "pushed clean chain with " << trackseed.size_cluster_keys() << " clusters" << std::endl;
    }
  }

  return clean_chains;
}

void PHCASeeding::publishSeeds(const std::vector<TrackSeed_v2>& seeds)
{
  for (const auto& seed : seeds)
  {
    auto pseed = std::make_unique<TrackSeed_v2>(seed);
    if (Verbosity() > 4)
    {
      pseed->identify();
    }
    _track_map->insert(pseed.get());
  }
}

int PHCASeeding::Setup(PHCompositeNode* topNode) // This is called by ::InitRun
{
  //  if(Verbosity()>0)
  std::cout << "Called Setup" << std::endl;
  if (Verbosity() > 0)
  {
    std::cout << "topNode:" << topNode << std::endl;
  }
  PHTrackSeeding::Setup(topNode);

  // geometry initialization
  int ret = InitializeGeometry(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  // tpc distortion correction
  m_dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (m_dcc)
  {
    std::cout << "PHCASeeding::Setup - found static TPC distortion correction container" << std::endl;
  }

  t_fill = std::make_unique<PHTimer>("t_fill");
  t_fill->stop();

  t_seed = std::make_unique<PHTimer>("t_seed");
  t_seed->stop();

  t_makebilinks = std::make_unique<PHTimer>("t_makebilinks");
  t_makebilinks->stop();

  t_makeseeds = std::make_unique<PHTimer>("t_makeseeds");
  t_makeseeds->stop();

  //  fcfg.set_rescale(1);
  std::unique_ptr<PHField> field_map;
  if (_use_const_field)
  {
    PHFieldConfigv2 fcfg(0, 0, _const_field);
    field_map = std::unique_ptr<PHField>(PHFieldUtility::BuildFieldMap(&fcfg));
  }
  else
  {
    PHFieldConfigv1 fcfg;
    fcfg.set_field_config(PHFieldConfig::FieldConfigTypes::Field3DCartesian);
    if (std::filesystem::path(m_magField).extension() != ".root")
    {
      m_magField = CDBInterface::instance()->getUrl(m_magField);
    }
    fcfg.set_filename(m_magField);
    field_map = std::unique_ptr<PHField>(PHFieldUtility::BuildFieldMap(&fcfg));
  }

  fitter = std::make_unique<ALICEKF>(topNode, _cluster_map, field_map.get(), _fieldDir, _min_clusters_per_track, _max_sin_phi, Verbosity());
  fitter->useConstBField(_use_const_field);
  fitter->setConstBField(_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->setFixedClusterError(0, _fixed_clus_err.at(0));
  fitter->setFixedClusterError(1, _fixed_clus_err.at(1));
  fitter->setFixedClusterError(2, _fixed_clus_err.at(2));
  
#if defined(_CLUSTER_LOG_TUPOUT_)
  std::cout << " Writing _CLUSTER_LOG_TUPOUT.root file " << std::endl;
  _f_clustering_process = new TFile("_CLUSTER_LOG_TUPOUT.root", "recreate");
  _tupclus_all         = new TNtuple("all",         "all clusters","event:layer:num:x:y:z");
  _tupclus_links       = new TNtuple("links",       "links","event:layer:updown01:x:y:z");
  _tupclus_bilinks     = new TNtuple("bilinks",     "bilinks","event:layer:topbot01:x:y:z");
  _tupclus_seeds       = new TNtuple("seeds",       "3 bilink seeds cores","event:layer:seed012:x:y:z");
  _tupclus_grown_seeds = new TNtuple("grown_seeds", "grown seeds", "event:layer:seednum05:x:y:z");
#endif


  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::End()
{
  if (Verbosity() > 0)
  {
    std::cout << "Called End " << std::endl;
  }

#if defined(_CLUSTER_LOG_TUPOUT_)
  _f_clustering_process->cd();
  _tupclus_all         ->Write();
  _tupclus_links       ->Write();
  _tupclus_bilinks     ->Write();
  _tupclus_seeds       ->Write();
  _tupclus_grown_seeds ->Write();
  _f_clustering_process->Close();
#endif

  return Fun4AllReturnCodes::EVENT_OK;
}
    // ---OLD CODE 1: SKIP_LAYERS---
  //  trackSeedKeyLists = tempSeedKeyLists;
  /*
    for(auto trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
    {
      bool reached_end = false;
      while(!reached_end)
      {
        TrkrDefs::cluskey trackHead = trackKeyChain->back();
        TrkrDefs::cluskey secondToLast = trackKeyChain->rbegin()[1];
        TrkrDefs::cluskey thirdToLast = trackKeyChain->rbegin()[2];
        auto& head_pos = globalPositions.at(trackHead);
        auto& sec_pos = globalPositions.at(secondToLast);
        auto& third_pos = globalPositions.at(thirdToLast);
        double dz_avg = ((head_pos.z()-sec_pos.z())+(sec_pos.z()-third_pos.z()))/2.;
        double dx1 = head_pos.x()-sec_pos.x();
        double dy1 = head_pos.y()-sec_pos.y();
        double dx2 = sec_pos.x()-third_pos.x();
        double dy2 = sec_pos.y()-third_pos.y();
        double ddx = dx1-dx2;
        double ddy = dy1-dy2;
        double new_dx = dx1+ddx;
        double new_dy = dy1+ddy;
        double new_x = head_pos.x()+new_dx;
        double new_y = head_pos.y()+new_dy;
        double new_z = head_pos.z()+dz_avg;
        std::cout << "(x,y,z) = (" << head_pos.x() << ", " << head_pos.y() << ", " << head_pos.z() << ")" << std::endl;
        unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - (_nlayers_intt + _nlayers_maps);
        std::cout << "layer " << trackHead_layer << std::endl;
        std::cout << "projected: (" << new_x << ", " << new_y << ", " << new_z << ")" << std::endl;
        TrkrDefs::cluskey nextCluster;
        double bestDist = 1e9;
        bool no_next_link = true;
        for(auto testlink = bilinks[trackHead_layer].begin(); testlink != bilinks[trackHead_layer].end(); ++testlink)
        {
          if((*testlink)[0].second==trackHead)
          {
            TrkrDefs::cluskey testCluster = (*testlink)[1].second;
            auto& test_pos = globalPositions.at(testCluster);
            std::cout << "test cluster: (" << test_pos.x() << ", " << test_pos.y() << ", " << test_pos.z() << ")" << std::endl;
            double distToNew = sqrt(square<double>(test_pos.x()-new_x)+square<double>(test_pos.y()-new_y)+square<double>(test_pos.z()-new_z));
            if(distToNew<bestDist)
            {
              std::cout << "current best" << std::endl;
              nextCluster = testCluster;
              bestDist = distToNew;
            }
            no_next_link = false;
          }
        }
        if(!no_next_link) trackKeyChain->push_back(nextCluster);
        if(no_next_link) reached_end = true;
      }
    }
  */
    // ---OLD CODE 0: SKIP_LAYERS---
    /*
    if(mode == skip_layers::on)
    {
      if(maxCosPlaneAngle > _cosTheta_limit)
      {
        // if no triplet is sufficiently linear, then it's likely that there's a missing cluster
        // repeat search but skip one layer below
        std::vector<pointKey> clustersTwoLayersBelow;
        QueryTree(_rtree,
                StartPhi-_neighbor_phi_width,
                StartEta-_neighbor_eta_width,
                (double) StartLayer - 2.5,
                StartPhi+_neighbor_phi_width,
                StartEta+_neighbor_eta_width,
                (double) StartLayer - 1.5,
                clustersTwoLayersBelow);
        std::vector<std::array<double,3>> delta_2below;
        delta_2below.clear();
        delta_2below.resize(clustersTwoLayersBelow.size());
        std::transform(clustersTwoLayersBelow.begin(),clustersTwoLayersBelow.end(),delta_2below.begin(),
          [&](pointKey BelowCandidate){
            const auto& belowpos = globalPositions.at(BelowCandidate.second);
            return std::array<double,3>{(belowpos(0))-StartX,
              (belowpos(1))-StartY,
              (belowpos(2))-StartZ};});
        for(size_t iAbove = 0; iAbove<delta_above.size(); ++iAbove)
        {
          for(size_t iBelow = 0; iBelow<delta_2below.size(); ++iBelow)
          {
            double dotProduct = delta_2below[iBelow][0]*delta_above[iAbove][0]+delta_2below[iBelow][1]*delta_above[iAbove][1]+delta_2below[iBelow][2]*delta_above[iAbove][2];
            double belowSqLength = sqrt(delta_2below[iBelow][0]*delta_2below[iBelow][0]+delta_2below[iBelow][1]*delta_2below[iBelow][1]+delta_2below[iBelow][2]*delta_2below[iBelow][2]);
            double aboveSqLength = sqrt(delta_above[iAbove][0]*delta_above[iAbove][0]+delta_above[iAbove][1]*delta_above[iAbove][1]+delta_above[iAbove][2]*delta_above[iAbove][2]);
            double cosPlaneAngle = dotProduct / (belowSqLength*aboveSqLength);
            if(cosPlaneAngle < maxCosPlaneAngle)
            {
              maxCosPlaneAngle = cosPlaneAngle;
              bestBelowCluster = fromPointKey(clustersTwoLayersBelow[iBelow]);
              bestAboveCluster = fromPointKey(ClustersAbove[iAbove]);
            }
          }
        }
        // if no triplet is STILL sufficiently linear, then do the same thing, but skip one layer above
        if(maxCosPlaneAngle > _cosTheta_limit)
        {
          std::vector<pointKey> clustersTwoLayersAbove;
          QueryTree(_rtree,
                  StartPhi-_neighbor_phi_width,
                  StartEta-_neighbor_eta_width,
                  (double) StartLayer + 1.5,
                  StartPhi+_neighbor_phi_width,
                  StartEta+_neighbor_eta_width,
                  (double) StartLayer + 2.5,
                  clustersTwoLayersAbove);
          std::vector<std::array<double,3>> delta_2above;
          delta_2above.clear();
          delta_2above.resize(clustersTwoLayersAbove.size());
          std::transform(clustersTwoLayersAbove.begin(),clustersTwoLayersAbove.end(),delta_2above.begin(),
            [&](pointKey AboveCandidate){
              const auto& abovepos = globalPositions.at(AboveCandidate.second);
              return std::array<double,3>{(abovepos(0))-StartX,
                (abovepos(1))-StartY,
                (abovepos(2))-StartZ};});
          for(size_t iAbove = 0; iAbove<delta_2above.size(); ++iAbove)
          {
            for(size_t iBelow = 0; iBelow<delta_below.size(); ++iBelow)
            {
              double dotProduct = delta_below[iBelow][0]*delta_2above[iAbove][0]+delta_below[iBelow][1]*delta_2above[iAbove][1]+delta_below[iBelow][2]*delta_2above[iAbove][2];
              double belowSqLength = sqrt(delta_below[iBelow][0]*delta_below[iBelow][0]+delta_below[iBelow][1]*delta_below[iBelow][1]+delta_below[iBelow][2]*delta_below[iBelow][2]);
              double aboveSqLength = sqrt(delta_2above[iAbove][0]*delta_2above[iAbove][0]+delta_2above[iAbove][1]*delta_2above[iAbove][1]+delta_2above[iAbove][2]*delta_2above[iAbove][2]);
              double cosPlaneAngle = dotProduct / (belowSqLength*aboveSqLength);
              if(cosPlaneAngle < maxCosPlaneAngle)
              {
                maxCosPlaneAngle = cosPlaneAngle;
                bestBelowCluster = fromPointKey(ClustersBelow[iBelow]);
                bestAboveCluster = fromPointKey(clustersTwoLayersAbove[iAbove]);
              }
            }
          }
        }
      }
    }
    */
