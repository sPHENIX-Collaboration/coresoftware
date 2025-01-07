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
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <phfield/PHFieldConfigv2.h>

#include <ffamodules/CDBInterface.h>

// trackbase_historic includes
#include <trackbase/ActsGeometry.h>
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
  if (Verbosity() > 2) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void) 0
#endif

#if defined(_PHCASEEDING_TIMER_OUT_)
#define PHCASEEDING_PRINT_TIME(timer, statement) \
  timer->stop();                                 \
  std::cout << " PHCASEEDING_PRINT_TIME: Time to " << statement << ": " << timer->elapsed() / 1000 << " s" << std::endl;
#else
#define PHCASEEDING_PRINT_TIME(timer, statement) (void) 0
#endif

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
  /* inline double get_eta(const Acts::Vector3& position) */
  /* { */
  /*   const double norm = std::sqrt(square(position.x()) + square(position.y()) + square(position.z())); */
  /*   return std::log((norm + position.z()) / (norm - position.z())) / 2; */
  /* } */

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
    float neighbor_z_width,
    float maxSinPhi)
  : PHTrackSeeding(name)
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _start_layer(start_layer)
  , _end_layer(end_layer)
  , _min_nhits_per_cluster(min_nhits_per_cluster)
  , _min_clusters_per_track(min_clusters_per_track)
  , _neighbor_phi_width(neighbor_phi_width)
  , _neighbor_z_width(neighbor_z_width)
  , _max_sin_phi(maxSinPhi)
{
}

int PHCASeeding::InitializeGeometry(PHCompositeNode* topNode)
{
  // geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No acts tracking geometry, can't proceed" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 PHCASeeding::getGlobalPosition(TrkrDefs::cluskey key, TrkrCluster* cluster) const
{
  return _pp_mode ? m_tGeometry->getGlobalPosition(key, cluster) : m_globalPositionWrapper.getGlobalPositionDistortionCorrected(key, cluster, 0);
}

void PHCASeeding::QueryTree(const bgi::rtree<PHCASeeding::pointKey, bgi::quadratic<16>>& rtree, double phimin, double z_min, double phimax, double z_max, std::vector<pointKey>& returned_values) const
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
    rtree.query(bgi::intersects(box(point(phimin, z_min), point(2 * M_PI, z_max))), std::back_inserter(returned_values));
    rtree.query(bgi::intersects(box(point(0., z_min), point(phimax, z_max))), std::back_inserter(returned_values));
  }
  else
  {
    rtree.query(bgi::intersects(box(point(phimin, z_min), point(phimax, z_max))), std::back_inserter(returned_values));
  }
}

std::pair<PHCASeeding::PositionMap, PHCASeeding::keyListPerLayer> PHCASeeding::FillGlobalPositions()
{
  keyListPerLayer ckeys;
  PositionMap cachedPositions;
  cachedPositions.reserve(_cluster_map->size());  // avoid resizing mid-execution

  for (const auto& hitsetkey : _cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto range = _cluster_map->getClusters(hitsetkey);
    for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
    {
      TrkrDefs::cluskey ckey = clusIter->first;
      TrkrCluster* cluster = clusIter->second;
      unsigned int layer = TrkrDefs::getLayer(ckey);

      if(cluster->getZSize()==1&&_reject_zsize1==true){
	continue;
      }
      if (layer < _start_layer || layer >= _end_layer)
      {
        if (Verbosity() > 2)
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

      ckeys[layer - _FIRST_LAYER_TPC].push_back(ckey);
      fill_tuple(_tupclus_all, 0, ckey, cachedPositions.at(ckey));
    }
  }
  return std::make_pair(cachedPositions, ckeys);
}

std::vector<PHCASeeding::coordKey> PHCASeeding::FillTree(bgi::rtree<PHCASeeding::pointKey, bgi::quadratic<16>>& _rtree, const PHCASeeding::keyList& ckeys, const PHCASeeding::PositionMap& globalPositions, const int layer)
{
  // Fill _rtree with the clusters in ckeys; remove duplicates, and return a vector of the coordKeys
  // Note that layer is only used for a cout statement
  int n_dupli = 0;
  std::vector<coordKey> coords;
  _rtree.clear();
  /* _rtree.reserve(ckeys.size()); */
  for (const auto& ckey : ckeys)
  {
    const auto& globalpos_d = globalPositions.at(ckey);
    const double clus_phi = get_phi(globalpos_d);
    const double clus_z = globalpos_d.z();
    if (Verbosity() > 5)
    {
      /* int layer = TrkrDefs::getLayer(ckey); */
      std::cout << "Found cluster " << ckey << " in layer " << layer << std::endl;
    }
    std::vector<pointKey> testduplicate;
    QueryTree(_rtree, clus_phi - 0.00001, clus_z - 0.00001, clus_phi + 0.00001, clus_z + 0.00001, testduplicate);
    if (!testduplicate.empty())
    {
      ++n_dupli;
      continue;
    }
    coords.push_back({{static_cast<float>(clus_phi), static_cast<float>(clus_z)}, ckey});
    t_fill->restart();
    _rtree.insert(std::make_pair(point(clus_phi, globalpos_d.z()), ckey));
    t_fill->stop();
  }
  if (Verbosity() > 5)
  {
    std::cout << "nhits in layer(" << layer << "): " << coords.size() << std::endl;
  }
  if (Verbosity() > 3)
  {
    std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
  }
  if (Verbosity() > 3)
  {
    std::cout << "number of duplicates : " << n_dupli << std::endl;
  }
  return coords;
}

int PHCASeeding::Process(PHCompositeNode* /*topNode*/)
{
  process_tupout_count();
  if (Verbosity() > 3)
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
  numberofseeds += FindSeedsWithMerger(globalPositions, ckeys);
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
  PHCASEEDING_PRINT_TIME(t_makebilinks, "init and make bilinks");

  t_makeseeds->restart();
  keyLists trackSeedKeyLists = FollowBiLinks(trackSeedPairs, bodyLinks, globalPositions);
  PHCASEEDING_PRINT_TIME(t_makeseeds, "make seeds");
  if (Verbosity() > 0)
  {
    t_makeseeds->stop();
    std::cout << "Time to make seeds: " << t_makeseeds->elapsed() / 1000 << " s" << std::endl;
  }
  std::vector<TrackSeed_v2> seeds = RemoveBadClusters(trackSeedKeyLists, globalPositions);

  publishSeeds(seeds);
  return seeds.size();
}

std::pair<PHCASeeding::keyLinks, PHCASeeding::keyLinkPerLayer> PHCASeeding::CreateBiLinks(const PHCASeeding::PositionMap& globalPositions, const PHCASeeding::keyListPerLayer& ckeys)
{
  keyLinks startLinks;        // bilinks at start of chains
  keyLinkPerLayer bodyLinks;  //  bilinks to build chains
                              //
  double cluster_find_time = 0;
  double rtree_query_time = 0;
  double transform_time = 0;
  double compute_best_angle_time = 0;
  double set_insert_time = 0;

  // there are three coord_array (only the current layer is used at a time,
  // but it is filled the same time as the _rtrees, which are used two at
  // a time -- the prior padplane row and the next padplain row
  std::array<std::vector<coordKey>, 3> coord_arr;
  std::array<std::unordered_set<keyLink>, 2> previous_downlinks_arr;
  std::array<std::unordered_set<TrkrDefs::cluskey>, 2> bottom_of_bilink_arr;

  // iterate from outer to inner layers
  const int inner_index = _start_layer - _FIRST_LAYER_TPC + 1;
  const int outer_index = _end_layer - _FIRST_LAYER_TPC - 2;

  // fill the current and prior row coord and ttrees for the first iteration
  int _index_above = (outer_index + 1) % 3;
  int _index_current = (outer_index) % 3;
  coord_arr[_index_above] = FillTree(_rtrees[_index_above], ckeys[outer_index + 1], globalPositions, outer_index + 1);
  coord_arr[_index_current] = FillTree(_rtrees[_index_current], ckeys[outer_index], globalPositions, outer_index);

  for (int layer_index = outer_index; layer_index >= inner_index; --layer_index)
  {
    // these lines of code will rotates through all three _rtree's in the array,
    // where the old lower becomes the new middle, the old middle the new upper,
    // and the old upper drops out and that _rtree is filled with the new lower
    const unsigned int LAYER = layer_index + _FIRST_LAYER_TPC;
    int index_above = (layer_index + 1) % 3;
    int index_current = (layer_index) % 3;
    int index_below = (layer_index - 1) % 3;

    coord_arr[index_below] = FillTree(_rtrees[index_below], ckeys[layer_index - 1], globalPositions, layer_index - 1);

    // NO DUPLICATES FOUND IN COORD_ARR

    auto& _rtree_above = _rtrees[index_above];
    const std::vector<coordKey>& coord = coord_arr[index_current];
    auto& _rtree_below = _rtrees[index_below];

    auto& curr_downlinks = previous_downlinks_arr[layer_index % 2];
    auto& last_downlinks = previous_downlinks_arr[(layer_index + 1) % 2];

    auto& curr_bottom_of_bilink = bottom_of_bilink_arr[layer_index % 2];
    auto& last_bottom_of_bilink = bottom_of_bilink_arr[(layer_index + 1) % 2];

    curr_downlinks.clear();
    curr_bottom_of_bilink.clear();

    // For all the clusters in coord, find nearest neighbors in the
    // above and below layers and make links
    // Any link to an above node which matches the same clusters
    // on the previous iteration (to a "below node") becomes a "bilink"
    // Check if this bilink links to a prior bilink or not

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
                StartPhi - dphi_per_layer[LAYER],
                StartZ - dZ_per_layer[LAYER],
                StartPhi + dphi_per_layer[LAYER],
                StartZ + dZ_per_layer[LAYER],
                ClustersBelow);

      FillTupWinLink(_rtree_below, StartCluster, globalPositions);

      QueryTree(_rtree_above,
                StartPhi - dphi_per_layer[LAYER + 1],
                StartZ - dZ_per_layer[LAYER + 1],
                StartPhi + dphi_per_layer[LAYER + 1],
                StartZ + dZ_per_layer[LAYER + 1],
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
      // calculate (delta_z_, delta_phi) vector for each neighboring cluster

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
      // (by maximizing the cos of the angle between the (delta_z_,delta_phi) vectors)
      // double minSumLengths = 1e9;
      std::unordered_set<TrkrDefs::cluskey> bestAboveClusters;
      for (size_t iAbove = 0; iAbove < delta_above.size(); ++iAbove)
      {
        for (size_t iBelow = 0; iBelow < delta_below.size(); ++iBelow)
        {
          // test for straightness of line just by taking the cos(angle) between the two vectors
          // use the sq as it is much faster than sqrt
          const auto& A = delta_below[iBelow];
          const auto& B = delta_above[iAbove];
          // calculate normalized dot product between two vectors
          const double A_len_sq = (A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
          const double B_len_sq = (B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
          const double dot_prod = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
          const double cos_angle_sq = dot_prod * dot_prod / A_len_sq / B_len_sq;  // also same as cos(angle), where angle is between two vectors
          FillTupWinCosAngle(ClustersAbove[iAbove].second, StartCluster.second, ClustersBelow[iBelow].second, globalPositions, cos_angle_sq, (dot_prod < 0.));

          constexpr double maxCosPlaneAngle = -0.95;
          constexpr double maxCosPlaneAngle_sq = maxCosPlaneAngle * maxCosPlaneAngle;
          if ((dot_prod < 0.) && (cos_angle_sq > maxCosPlaneAngle_sq))
          {
            // maxCosPlaneAngle = cos(angle);
            // minSumLengths = belowLength+aboveLength;
            curr_downlinks.insert({StartCluster.second, ClustersBelow[iBelow].second});
            bestAboveClusters.insert(ClustersAbove[iAbove].second);

            // fill the tuples for plotting
            fill_tuple(_tupclus_links, 0, StartCluster.second, globalPositions.at(StartCluster.second));
            fill_tuple(_tupclus_links, -1, ClustersBelow[iBelow].second, globalPositions.at(ClustersBelow[iBelow].second));
            fill_tuple(_tupclus_links, 1, ClustersAbove[iAbove].second, globalPositions.at(ClustersAbove[iAbove].second));
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

      for (auto cluster : bestAboveClusters)
      {
        keyLink uplink = std::make_pair(cluster, StartCluster.second);

        if (last_downlinks.find(uplink) != last_downlinks.end())
        {
          // this is a bilink
          const auto& key_top = uplink.first;
          const auto& key_bot = uplink.second;
          curr_bottom_of_bilink.insert(key_bot);
          fill_tuple(_tupclus_bilinks, 0, key_top, globalPositions.at(key_top));
          fill_tuple(_tupclus_bilinks, 1, key_bot, globalPositions.at(key_bot));

          if (last_bottom_of_bilink.find(key_top) == last_bottom_of_bilink.end())
          {
            startLinks.push_back(std::make_pair(key_top, key_bot));
          }
          else
          {
            bodyLinks[layer_index + 1].push_back(std::make_pair(key_top, key_bot));
          }
        }
      }  // end loop over all up-links
    }    // end loop over start clusters

    t_seed->stop();
    set_insert_time += t_seed->elapsed();
    t_seed->restart();
    LogDebug(" max collinearity: " << maxCosPlaneAngle << std::endl);
  }  // end loop over layers (to make links)

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

PHCASeeding::keyLists PHCASeeding::FollowBiLinks(const PHCASeeding::keyLinks& trackSeedPairs, const PHCASeeding::keyLinkPerLayer& bilinks, const PHCASeeding::PositionMap& globalPositions) const
{
  // form all possible starting 3-cluster tracks (we need that to calculate curvature)
  keyLists seeds;
  for (auto& startLink : trackSeedPairs)
  {
    TrkrDefs::cluskey trackHead = startLink.second;
    unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - _FIRST_LAYER_TPC;
    // the following call with get iterators to all bilinks which match the head
    for (const auto& matchlink : bilinks[trackHead_layer])
    {
      /* auto matched_links = std::equal_range(bilinks[trackHead_layer].begin(), bilinks[trackHead_layer].end(), trackHead, CompKeyToBilink()); */
      /* for (auto matchlink = matched_links.first; matchlink != matched_links.second; ++matchlink) */
      /* { */
      if (matchlink.first != trackHead)
      {
        continue;
      }
      keyList trackSeedTriplet;
      trackSeedTriplet.push_back(startLink.first);
      trackSeedTriplet.push_back(startLink.second);
      trackSeedTriplet.push_back(matchlink.second);
      seeds.push_back(trackSeedTriplet);

      fill_tuple(_tupclus_seeds, 0, startLink.first, globalPositions.at(startLink.first));
      fill_tuple(_tupclus_seeds, 1, startLink.second, globalPositions.at(startLink.second));
      fill_tuple(_tupclus_seeds, 2, matchlink.second, globalPositions.at(matchlink.second));
    }
  }

  // - grow every seed in the seedlist, up to the maximum number of clusters per seed
  // - the algorithm is that every cluster is allowed to be used by any number of chains, so there is no penalty in which order they are added

  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "starting cluster finding time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  t_seed->restart();
  // assemble track cluster chains from starting cluster keys (ordered from outside in)

  // std::cout << "STARTING SEED ASSEMBLY" << std::endl;

  // The algorithm is to add keylinks to every chain of bilinks (seed) that wants them
  // It is not first come first serve
  // Therefore, follow each chain to the end
  // If there are possible multiple links to add to a single chain, optionally split the chain to
  // follow all possibile links (depending on the input parameter _split_seeds)

  //  std::vector<keyList> seeds;
  // std::vector<keyList> tempSeedKeyLists = seeds;
  // seeds.clear();
  if (seeds.size() == 0)
  {
    return seeds;
  }

  int nsplit_chains = -1;

  // positions of the seed being following
  std::array<float, 4> phi{}, R{}, Z{};
  keyLists grown_seeds;

  while (seeds.size() > 0)
  {
    keyLists split_seeds{};  // to collect when using split tracks
    for (auto& seed : seeds)
    {
      // grow the seed to the maximum length allowed
      bool first_link = true;
      bool done_growing = (seed.size() >= _max_clusters_per_seed);
      keyList head_keys = {seed.back()};
      /* keyList head_keys = { seed.back() }; // heads of the seed */

      while (!done_growing)
      {
        // Get all bilinks which fit to the head of the chain
        unsigned int iL = TrkrDefs::getLayer(head_keys[0]) - _FIRST_LAYER_TPC;
        keySet link_matches{};
        for (const auto& head_key : head_keys)
        {
          // also possible to sort the links and use a sorted search like:
          // auto matched_links = std::equal_range(bilinks[trackHead_layer].begin(), bilinks[trackHead_layer].end(), trackHead, CompKeyToBilink());
          // for (auto link = matched_links.first; link != matched_links.second; ++link)
          for (auto& link : bilinks[iL])
          {  // iL for "Index of Layer"
            if (link.first == head_key)
            {
              link_matches.insert(link.second);
            }
          }
        }

        // find which link_matches pass the dZdR and d2phidr2 cuts
        keyList passing_links{};
        for (const auto link : link_matches)
        {  // iL for "Index of Layer"
          // see if the link passes the growth cuts
          if (first_link)
          {
            first_link = false;
            for (int i = 1; i < 4; ++i)
            {
              const auto& pos = globalPositions.at(seed.rbegin()[i - 1]);
              const auto x = pos.x();
              const auto y = pos.y();
              int index = (iL + i) % 4;
              Z[index] = pos.z();
              phi[index] = atan2(y, x);
              R[index] = sqrt(x * x + y * y);
            }
          }

          // get the data for the new link
          const auto& pos = globalPositions.at(link);
          const auto x = pos.x();
          const auto y = pos.y();
          const auto z = pos.z();

          const int i0 = (iL + 0) % 4;
          const int i1 = (iL + 1) % 4;
          const int i2 = (iL + 2) % 4;
          const int i3 = (iL + 3) % 4;

          phi[i0] = atan2(y, x);
          R[i0] = sqrt(x * x + y * y);
          Z[i0] = z;

          // see if it is possible matching link
          if (_split_seeds)
          {
            FillTupWinGrowSeed(seed, {head_keys[0], link}, globalPositions);
          }
          const float dZ_12 = Z[i1] - Z[i2];
          const float dZ_01 = Z[i0] - Z[i1];
          const float dR_12 = R[i1] - R[i2];
          const float dR_01 = R[i0] - R[i1];
          const float dZdR_01 = dZ_01 / dR_01;
          const float dZdR_12 = dZ_12 / dR_12;

          if (fabs(dZdR_01 - dZdR_12) > _clusadd_delta_dzdr_window)
          {
            continue;
          }
          const float dphi_01 = phi[i0] - phi[i1];
          const float dphi_12 = phi[i1] - phi[i2];
          const float dphi_23 = phi[i2] - phi[i3];
          const float dR_23 = R[i2] - R[i3];
          const float d2phidr2_01 = dphi_01 / dR_01 / dR_01 - dphi_12 / dR_12 / dR_12;
          const float d2phidr2_12 = dphi_12 / dR_12 / dR_12 - dphi_23 / dR_23 / dR_23;
          if (fabs(d2phidr2_01 - d2phidr2_12) > _clusadd_delta_dphidr2_window)
          {
            continue;
          }
          passing_links.push_back(link);
        }  // end loop over all bilinks in new layer

        if (_split_seeds)
        {
          fill_split_chains(seed, passing_links, globalPositions, nsplit_chains);
        }

        // grow the chain appropriately
        switch (passing_links.size())
        {
        case 0:
          done_growing = true;
          break;
        case 1:
          seed.push_back(passing_links[0]);
          if (seed.size() >= _max_clusters_per_seed)
          {
            done_growing = true;
          }  // this seed is done growing
          head_keys = {passing_links[0]};
          break;
        default:  // more than one matched cluster
          if (_split_seeds)
          {
            // there are multiple matching clusters
            // if we are splitting seeds, then just push back each of the matched
            // to the back of the seeds to grow on their own
            for (unsigned int i = 1; i < passing_links.size(); ++i)
            {
              keyList newseed = {seed.begin(), seed.end()};
              newseed.push_back(passing_links[i]);
              split_seeds.push_back(newseed);
            }
            seed.push_back(passing_links[0]);
            if (seed.size() >= _max_clusters_per_seed)
            {
              done_growing = true;
            }
            head_keys = {passing_links[0]};
          }
          else
          {
            // multiple seeds matched. get the average position to put into
            // Z, phi, and R (of [iL]), and pass all the links to find the next cluster
            float avg_x = 0;
            float avg_y = 0;
            float avg_z = 0;
            for (const auto& link : passing_links)
            {
              const auto& pos = globalPositions.at(link);
              avg_x += pos.x();
              avg_y += pos.y();
              avg_z += pos.z();
            }
            avg_x /= passing_links.size();
            avg_y /= passing_links.size();
            avg_z /= passing_links.size();
            phi[iL % 4] = atan2(avg_y, avg_x);
            R[iL % 4] = sqrt(avg_x * avg_x + avg_y * avg_y);
            Z[iL % 4] = avg_z;
            head_keys = passing_links;  // will try and grow from this position
          }                             // end of logic for processing passing seeds
          break;
        }  // end of seed length switch
      }    // end of seed growing loop: if (!done_growing)
      if (seed.size() >= _min_clusters_per_seed)
      {
        grown_seeds.push_back(seed);
        fill_tuple_with_seed(_tupclus_grown_seeds, seed, globalPositions);
      }
    }  // end of loop over seeds
    seeds.clear();
    for (const auto& seed : split_seeds)
    {
      seeds.push_back(seed);
    }
    /* seeds = split_seeds; */
  }  // end of looping over all seeds

  // old code block move to end of code under the title: "---OLD CODE 1: SKIP_LAYERS---"
  t_seed->stop();
  if (Verbosity() > 1)
  {
    std::cout << "keychain assembly time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  t_seed->restart();
  LogDebug(" track key chains assembled: " << trackSeedKeyLists.size() << std::endl);
  LogDebug(" track key chain lengths: " << std::endl);
  return grown_seeds;
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
    if (Verbosity() > 3)
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
    for (const auto& key : chain)
    {
      trackseed.insert_cluster_key(key);
    }
    clean_chains.push_back(trackseed);
    if (Verbosity() > 2)
    {
      std::cout << "pushed clean chain with " << trackseed.size_cluster_keys() << " clusters" << std::endl;
    }
  }

  return clean_chains;
}

void PHCASeeding::publishSeeds(const std::vector<TrackSeed_v2>& seeds) const
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

int PHCASeeding::Setup(PHCompositeNode* topNode)  // This is called by ::InitRun
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

  // timing
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
  fitter->setNeonFraction(Ne_frac);
  fitter->setArgonFraction(Ar_frac);
  fitter->setCF4Fraction(CF4_frac);
  fitter->setNitrogenFraction(N2_frac);
  fitter->setIsobutaneFraction(isobutane_frac);
  fitter->useConstBField(_use_const_field);
  fitter->setConstBField(_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->setFixedClusterError(0, _fixed_clus_err.at(0));
  fitter->setFixedClusterError(1, _fixed_clus_err.at(1));
  fitter->setFixedClusterError(2, _fixed_clus_err.at(2));

  PHG4TpcCylinderGeomContainer* geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cerr << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  for (int i = 8; i <= 54; ++i)
  {
    const float rad_0 = geom_container->GetLayerCellGeom(i - 1)->get_radius();
    const float rad_1 = geom_container->GetLayerCellGeom(i)->get_radius();
    const float delta_rad = rad_1 - rad_0;

    // sector boundaries between 22-23 and 39-40

    dZ_per_layer[i] = _neighbor_z_width * delta_rad;
    dphi_per_layer[i] = _neighbor_phi_width * delta_rad;
  }

#if defined(_PHCASEEDING_CLUSTERLOG_TUPOUT_)
  std::cout << " Writing _CLUSTER_LOG_TUPOUT.root file " << std::endl;
  _f_clustering_process = new TFile("_CLUSTER_LOG_TUPOUT.root", "recreate");
  _tupclus_all = new TNtuple("all", "all clusters", "event:layer:num:x:y:z");
  _tupclus_links = new TNtuple("links", "links", "event:layer:updown01:x:y:z:delta_z:delta_phi");
  _tupclus_bilinks = new TNtuple("bilinks", "bilinks", "event:layer:topbot01:x:y:z");
  _tupclus_seeds = new TNtuple("seeds", "3 bilink seeds cores", "event:layer:seed012:x:y:z");
  _tupclus_grown_seeds = new TNtuple("grown_seeds", "grown seeds", "event:layer:seednum05:x:y:z");
  _tupwin_link = new TNtuple("win_link", "neighbor clusters considered to make links", "event:layer0:x0:y0:z0:layer1:x1:y1:z1:dphi:dz");
  _tupwin_cos_angle = new TNtuple("win_cos_angle", "cos angle to make links", "event:layer0:x0:y0:z0:layer1:x1:y1:z1:layer2:x2:y2:z2:cos_angle");
  _tupwin_seed23 = new TNtuple("win_seed23", "xyL for points 1 and 2", "event:layer2:x2:y2:z2:layer3:x3:y3:z3");
  _tupwin_seedL1 = new TNtuple("win_seedL1", "xyL+link stats for points 0 and Link",
                               "event:layerL:xL:yL:zL:layer1:x1:y1:z1:dzdr_12:dzdr_L1:delta_dzdr_12_L1:d2phidr2_123:d2phidr2_L12:delta_d2phidr2");

  _search_windows = new TNtuple("search_windows", "windows used in algorithm to seed clusters",
                                "DelZ_ClSearch:DelPhi_ClSearch:start_layer:end_layer:dzdr_ClAdd:dphidr2_ClAdd");

  _tup_chainfork = new TNtuple("chainfork", "chain with multiple links, which if forking", "event:nchain:layer:x:y:z:dzdr:d2phidr2:nlink:nlinks");  // nlinks to add, 0 ... nlinks
  _tup_chainbody = new TNtuple("chainbody", "chain body with multiple link options", "event:nchain:layer:x:y:z:dzdr:d2phidr2:nlink:nlinks");        // nlinks in chain being added to will be 0, 1, 2 ... working backward from the fork -- dZ and dphi are dropped for final links as necessary
#endif

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::End()
{
  if (Verbosity() > 0)
  {
    std::cout << "Called End " << std::endl;
  }
  write_tuples();  // if defined _PHCASEEDING_CLUSTERLOG_TUPOUT_
  return Fun4AllReturnCodes::EVENT_OK;
}

#if defined(_PHCASEEDING_CHAIN_FORKS_)

void PHCASeeding::fill_split_chains(const PHCASeeding::keyList& seed, const PHCASeeding::keyList& add_links, const PHCASeeding::PositionMap& pos, int& n_tupchains) const
{
  if (add_links.size() < 2) return;
  n_tupchains += 1;

  // fill the chain leading to the forks
  int nlinks = seed.size();

  bool has_0 = false;
  float x0{0.}, y0{0.}, z0{0.}, r0{0.}, phi0{0.};

  bool has_1 = false;
  float /*x1, y1,*/ z1{0.}, r1{0.}, phi1{0.};

  bool has_2 = false;
  float /*x2, y2, z2,*/ r2{0.}, phi2{0.};

  int index = seed.size();  // index will count backwards from the end of the chain
                            // dR is not calculated for the final (outermost layer) link
                            // dPhi is not calculated for the final two (outmost layer) link

  float dphidr01 = -1000.;
  for (const auto& link : seed)
  {
    index -= 1;
    if (has_1)
    {
      has_2 = true;
      /*x2=x1; y2=y1; z2=z1;*/ r2 = r1;
      phi2 = phi1;
    }
    if (has_0)
    {
      has_1 = true;
      /*x1=x0; y1=y0;*/ z1 = z0;
      r1 = r0;
      phi1 = phi0;
    }
    auto link_pos = pos.at(link);
    has_0 = true;
    x0 = link_pos.x();
    y0 = link_pos.y();
    z0 = link_pos.z();
    phi0 = atan2(y0, x0);
    r0 = sqrt(x0 * x0 + y0 * y0);

    if (!has_1)
    {
      _tup_chainbody->Fill(_tupout_count, n_tupchains, TrkrDefs::getLayer(link), x0, y0, z0, -1000., -1000., index, nlinks);
      continue;
    }
    float dzdr = (z0 - z1) / (r0 - r1);
    if (!has_2)
    {
      _tup_chainbody->Fill(_tupout_count, n_tupchains, TrkrDefs::getLayer(link), x0, y0, z0, dzdr, -1000., index, nlinks);
      continue;
    }
    float dphi01 = std::fmod(phi1 - phi0, M_PI);
    float dphi12 = std::fmod(phi2 - phi1, M_PI);
    float dr_01 = r1 - r0;
    float dr_12 = r2 - r1;
    dphidr01 = dphi01 / dr_01 / dr_01;
    float d2phidr2 = dphidr01 - dphi12 / dr_12 / dr_12;
    _tup_chainbody->Fill(_tupout_count, n_tupchains, TrkrDefs::getLayer(link), x0, y0, z0, dzdr, d2phidr2, index, nlinks);
  }

  // now fill a chain of the possible added seeds
  index = -1;
  nlinks = add_links.size();

  for (const auto& link : add_links)
  {
    index += 1;
    auto link_pos = pos.at(link);
    float xt = link_pos.x();
    float yt = link_pos.y();
    float zt = link_pos.z();
    float phit = atan2(yt, xt);
    float rt = sqrt(xt * xt + yt * yt);
    float dr_t0 = rt - r0;

    float dzdr = (zt - z0) / (rt - r0);
    float dphit0 = std::fmod(phit - phi0, M_PI);

    float d2phidr2 = dphit0 / dr_t0 / dr_t0 - dphidr01;
    _tup_chainfork->Fill(_tupout_count, n_tupchains, TrkrDefs::getLayer(link), xt, yt, zt, dzdr, d2phidr2, index, nlinks);
  }
}

#else
void PHCASeeding::fill_split_chains(const PHCASeeding::keyList& /*chain*/, const PHCASeeding::keyList& /*links*/, const PHCASeeding::PositionMap& /*pos*/, int& /*nchains*/) const {};
#endif

#if defined(_PHCASEEDING_CLUSTERLOG_TUPOUT_)
void PHCASeeding::write_tuples()
{
  _f_clustering_process->cd();
  _tupclus_all->Write();
  _tupclus_links->Write();
  _tupclus_bilinks->Write();
  _tupclus_seeds->Write();
  _tupclus_grown_seeds->Write();
  _tupwin_link->Write();
  _tupwin_cos_angle->Write();
  _tupwin_seed23->Write();
  _tupwin_seedL1->Write();
  _search_windows->Write();
  _tup_chainbody->Write();
  _tup_chainfork->Write();
  _f_clustering_process->Close();
}

void PHCASeeding::fill_tuple(TNtuple* tup, float val, TrkrDefs::cluskey key, const Acts::Vector3& pos) const
{
  tup->Fill(_tupout_count, TrkrDefs::getLayer(key), val, pos[0], pos[1], pos[2]);
}

void PHCASeeding::fill_tuple_with_seed(TNtuple* tup, const PHCASeeding::keyList& seed, const PHCASeeding::PositionMap& pos) const
{
  for (unsigned int i = 0; i < seed.size(); ++i)
  {
    fill_tuple(tup, (float) i, seed[i], pos.at(seed[i]));
  }
}

void PHCASeeding::process_tupout_count()
{
  _tupout_count += 1;
  if (_tupout_count != 0) return;
  _search_windows->Fill(_neighbor_z_width, _neighbor_phi_width, _start_layer, _end_layer, _clusadd_delta_dzdr_window, _clusadd_delta_dphidr2_window);
}

void PHCASeeding::FillTupWinLink(bgi::rtree<PHCASeeding::pointKey, bgi::quadratic<16>>& _rtree_below, const PHCASeeding::coordKey& StartCluster, const PHCASeeding::PositionMap& globalPositions) const
{
  double StartPhi = StartCluster.first[0];
  const auto& P0 = globalPositions.at(StartCluster.second);
  double StartZ = P0(2);
  // Fill TNTuple _tupwin_link
  std::vector<pointKey> ClustersBelow;
  QueryTree(_rtree_below,
            StartPhi - 1.,
            StartZ - 20.,
            StartPhi + 1.,
            StartZ + 20.,
            ClustersBelow);

  for (const auto& pkey : ClustersBelow)
  {
    const auto P1 = globalPositions.at(pkey.second);
    double dphi = bg::get<0>(pkey.first) - StartPhi;
    double dZ = P1(2) - StartZ;
    _tupwin_link->Fill(_tupout_count, TrkrDefs::getLayer(StartCluster.second), P0(0), P0(1), P0(2), TrkrDefs::getLayer(pkey.second), P1(0), P1(1), P1(2), dphi, dZ);
  }
}

void PHCASeeding::FillTupWinCosAngle(const TrkrDefs::cluskey A, const TrkrDefs::cluskey B, const TrkrDefs::cluskey C, const PHCASeeding::PositionMap& globalPositions, double cos_angle_sq, bool isneg) const
{
  // A is top cluster, B the middle, C the bottom
  // a,b,c are the positions

  auto a = globalPositions.at(A);
  auto b = globalPositions.at(B);
  auto c = globalPositions.at(C);

  _tupwin_cos_angle->Fill(_tupout_count,
                          TrkrDefs::getLayer(A), a[0], a[1], a[2],
                          TrkrDefs::getLayer(B), b[0], b[1], b[2],
                          TrkrDefs::getLayer(C), c[0], c[1], c[2],
                          (isneg ? -1 : 1) * sqrt(cos_angle_sq));
}

void PHCASeeding::FillTupWinGrowSeed(const PHCASeeding::keyList& seed, const PHCASeeding::keyLink& link, const PHCASeeding::PositionMap& globalPositions) const
{
  TrkrDefs::cluskey trackHead = seed.back();
  auto& head_pos = globalPositions.at(trackHead);
  auto& prev_pos = globalPositions.at(seed.rbegin()[1]);
  float x1 = head_pos.x();
  float y1 = head_pos.y();
  float z1 = head_pos.z();
  float x2 = prev_pos.x();
  float y2 = prev_pos.y();
  float z2 = prev_pos.z();
  float dr_12 = sqrt(x1 * x1 + y1 * y1) - sqrt(x2 * x2 + y2 * y2);
  /* TrkrDefs::cluskey testCluster = link.second; */
  auto& test_pos = globalPositions.at(link.second);
  float xt = test_pos.x();
  float yt = test_pos.y();
  float zt = test_pos.z();
  float dr_t1 = sqrt(xt * xt + yt * yt) - sqrt(x1 * x1 + y1 * y1);
  float dzdr_12 = (z1 - z2) / dr_12;
  float dzdr_t1 = (zt - z1) / dr_t1;
  // if (fabs(dzdr_12 - dzdr_t1) > _clusadd_delta_dzdr_window)) // then fail this link

  auto& third_pos = globalPositions.at(seed.rbegin()[2]);
  float x3 = third_pos.x();
  float y3 = third_pos.y();
  float z3 = third_pos.z();
  float dr_23 = sqrt(x2 * x2 + y2 * y2) - sqrt(x3 * x3 + y3 * y3);
  float phi1 = atan2(y1, x1);
  float phi2 = atan2(y2, x2);
  float phi3 = atan2(y3, x3);
  float dphi12 = std::fmod(phi1 - phi2, M_PI);
  float dphi23 = std::fmod(phi2 - phi3, M_PI);
  float d2phidr2_123 = dphi12 / (dr_12 * dr_12) - dphi23 / (dr_23 * dr_23);
  float dphit1 = std::fmod(atan2(yt, xt) - atan2(y1, x1), M_PI);
  float d2phidr2_t12 = dphit1 / (dr_t1 * dr_t1) - dphi12 / (dr_12 * dr_12);
  _tupwin_seed23->Fill(_tupout_count,
                       (TrkrDefs::getLayer(seed.rbegin()[1])), x2, y2, z2,
                       (TrkrDefs::getLayer(seed.rbegin()[2])), x3, y3, z3);
  _tupwin_seedL1->Fill(_tupout_count,
                       (TrkrDefs::getLayer(link.second)), xt, yt, zt,
                       (TrkrDefs::getLayer(seed.back())), x1, y1, z1,
                       dzdr_12, dzdr_t1, fabs(dzdr_12 - dzdr_t1),
                       d2phidr2_123, d2phidr2_t12, fabs(d2phidr2_123 - d2phidr2_t12));
}
#else
void PHCASeeding::write_tuples(){};
void PHCASeeding::fill_tuple(TNtuple* /**/, float /**/, TrkrDefs::cluskey /**/, const Acts::Vector3& /**/) const {};
void PHCASeeding::fill_tuple_with_seed(TNtuple* /**/, const PHCASeeding::keyList& /**/, const PHCASeeding::PositionMap& /**/) const {};
void PHCASeeding::process_tupout_count(){};
void PHCASeeding::FillTupWinLink(bgi::rtree<PHCASeeding::pointKey, bgi::quadratic<16>>& /**/, const PHCASeeding::coordKey& /**/, const PHCASeeding::PositionMap& /**/) const {};
void PHCASeeding::FillTupWinCosAngle(const TrkrDefs::cluskey /**/, const TrkrDefs::cluskey /**/, const TrkrDefs::cluskey /**/, const PHCASeeding::PositionMap& /**/, double /**/, bool /**/) const {};
void PHCASeeding::FillTupWinGrowSeed(const PHCASeeding::keyList& /**/, const PHCASeeding::keyLink& /**/, const PHCASeeding::PositionMap& /**/) const {};
#endif  // defined _PHCASEEDING_CLUSTERLOG_TUPOUT_

// ---OLD CODE 1: SKIP_LAYERS---
//  trackSeedKeyLists = tempSeedKeyLists;
/*
  for(auto trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)

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
