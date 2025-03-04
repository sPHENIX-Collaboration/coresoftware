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

  inline float wrap_dphi(float a, float b) {
    float   _dphi = b-a;
    return (_dphi < -M_PI) ? _dphi += 2*M_PI
         : (_dphi >  M_PI) ? _dphi -= 2*M_PI
         : _dphi;
  }

  inline std::array<double,3> keyPointDiff(const PHCASeeding::keyPoint& a, const PHCASeeding::keyPoint& b) {
    return {b.x-a.x, b.y-a.y, b.z-a.z};
  }

  void print_reset_time(std::unique_ptr<PHTimer>& timer, int verbosity, int threshhold, const std::string& msg) {
    timer->stop();
    if (verbosity>=threshhold) {
      std::cout << " || " << timer->elapsed()/1000. << " s || " << msg << std::endl;
    }
    timer->restart();
  }

}  // namespace


PHCASeeding::Checker_dphidz::Checker_dphidz(
  const float& _clusadd_delta_z_window,
  const float& _clusadd_delta_phi_window,
  keyPtrList& seed)
: delta_dzdr_window {_clusadd_delta_z_window},
  delta_dphidr2_window {_clusadd_delta_phi_window}
{
  // add the last three clusters from the seed triplet
  for (int i=0;i<3;++i) {
    const auto&  p = *(seed.end()-3+i);
    const float x = p->x;
    const float y = p->y;
    z[i] = p->z;
    phi[i] = p->phi;
    R[i] = sqrt(x*x+y*y);
    if (i>0) {
      dR[i] = R[i]-R[i-1];
      dZdR[i] = (z[i]-z[i-1])/dR[i];
      auto dphi = wrap_dphi(phi[i-1],phi[i]);
      dphidR2[i] = dphi/dR[i]/dR[i];
    }
  }
}

bool PHCASeeding::Checker_dphidz::check_cluster(PHCASeeding::keyPtr p)
{
  update(p);
  return 
       (fabs(dZdR[i3]-dZdR[i2]) <= delta_dzdr_window)
   &&  (fabs(dphidR2[i3]-2*dphidR2[i2]+dphidR2[i1]) <= delta_dphidr2_window);
}

void PHCASeeding::Checker_dphidz::update(const PHCASeeding::keyPtr p) {
  const float x = p->x;
  const float y = p->y;
  z[i3] = p->z;
  phi[i3] = p->phi;
  R[i3] = sqrt(x*x+y*y);
  dR[i3] = R[i3]-R[i2];
  dZdR[i3] = (z[i3]-z[i2])/dR[i3];
  auto dphi = wrap_dphi(phi[i2],phi[i3]);
  dphidR2[i3] = dphi/dR[i3]/dR[i3];
}

void PHCASeeding::Checker_dphidz::update(const PHCASeeding::keyPtrList& ptrs) {
  // get the average values of z, phi, R, and update with those
  double _x = std::accumulate(ptrs.begin(), ptrs.end(), 0., 
    [](double sum, const keyPtr& p) { return sum + p->x; });
  double _y = std::accumulate(ptrs.begin(), ptrs.end(), 0., 
    [](double sum, const keyPtr& p) { return sum + p->y; });
  double _z = std::accumulate(ptrs.begin(), ptrs.end(), 0., 
    [](double sum, const keyPtr& p) { return sum + p->z; });
  double n = ptrs.size();
  _x /= n;
  _y /= n;
  _z /= n;
  z[i3] = _z;
  phi[i3] = atan2(_y,_x);
  R[i3] = sqrt(_x*_x+_y*_y);
  dR[i3] = R[i3]-R[i2];
  dZdR[i3] = (z[i3]-z[i2])/dR[i3];
  dphidR2[i3] = (phi[i3]-phi[i2])/dR[i3]/dR[i3];
}

void PHCASeeding::Checker_dphidz::add_cluster(const PHCASeeding::keyPtr p) {
  if (p) { update(p); }
  ++index;
  i1 = (index+1)%4;
  i2 = (index+2)%4;
  i3 = (index+3)%4;
}

void PHCASeeding::Checker_dphidz::add_clusters(const keyPtrList& ptrs) {
  update(ptrs);
  ++index;
  i1 = (index+1)%4;
  i2 = (index+2)%4;
  i3 = (index+3)%4;
}

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

void PHCASeeding::QueryTree(
    const PHCASeeding::boost_rtree& rtree, 
    const PHCASeeding::keyPtr& p, 
    const double& dphi, const double& dz,  PHCASeeding::rtreePairList& returned_values) const
{
  double phimin = p->phi-dphi;
  double phimax = p->phi+dphi;
  double zmin = p->z-dz;
  double zmax = p->z+dz;
  bool query_both_ends = false;
  if (phimin < -M_PI)
  {
    query_both_ends = true;
    phimin += 2 * M_PI;
  }
  if (phimax > M_PI)
  {
    query_both_ends = true;
    phimax -= 2 * M_PI;
  }
  if (query_both_ends)
  {
    rtree.query(bgi::intersects(box(point(phimin, zmin), point(M_PI, zmax))), std::back_inserter(returned_values));
    rtree.query(bgi::intersects(box(point(-M_PI, zmin), point(phimax, zmax))), std::back_inserter(returned_values));
  }
  else
  {
    rtree.query(bgi::intersects(box(point(phimin, zmin), point(phimax, zmax))), std::back_inserter(returned_values));
  }
}

PHCASeeding::keyPtrArr PHCASeeding::GetKeyPoints() {
  _cluster_pts.clear();
  _cluster_pts.reserve(_cluster_map->size());

  // Fill the vector entirely for all data_points
  // make accessor to iterate over point in each layer
  std::array<int, _NLAYERS_TPC> cnt_per_layer{};

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
      _cluster_pts.emplace_back(ckey, globalpos_d);
      cnt_per_layer[layer-_FIRST_LAYER_TPC]++;
      fill_tuple(_tupclus_all, 0, &_cluster_pts.back());
    }
  }

  keyPtrArr ptsPerLayer{};
  for (int i=0;i<_NLAYERS_TPC;++i) {
    ptsPerLayer[i].reserve(cnt_per_layer[i]);
  }
  for (auto& p : _cluster_pts) {
    ptsPerLayer[TrkrDefs::getLayer(p.key)-_FIRST_LAYER_TPC].emplace_back(&p);
  }
  return ptsPerLayer;
}


PHCASeeding::keyPtrList PHCASeeding::FillTree(PHCASeeding::boost_rtree& _rtree, const PHCASeeding::keyPtrList& ptrs)
{
  // Fill _rtree with locations in the pointKeys; skip duplicates; and return a vector pointKeyIterators
  // Note: layer used only for cout statement
  int n_dupli = 0;
  keyPtrList iters_out; // vector<keyPointIer>
  iters_out.reserve(ptrs.size());
  _rtree.clear();
  /* _rtree.reserve(ckeys.size()); */
  t_fill->restart();
  for (const auto& p : ptrs)
  {
    if (Verbosity() > 5)
    {
      int layer = TrkrDefs::getLayer(p->key);
      std::cout << "Found cluster " << p->key << " in layer " << layer << std::endl;
    }
    rtreePairList testduplicate;
    QueryTree(_rtree, p, 0.00001, 0.00001, testduplicate);
    if (!testduplicate.empty())
    {
      ++n_dupli;
      continue;
    }
    iters_out.push_back(p);
    _rtree.insert(std::make_pair(point(p->phi, p->z), p));
  }
  t_fill->stop();
  if (Verbosity() > 5)
  {
    int layer = -1;
    if (ptrs.size()>0) { TrkrDefs::getLayer(ptrs[0]->key); };
    std::cout << "nhits in layer(" << layer << "): " << iters_out.size() << std::endl;
  }
  if (Verbosity() > 3)
  {
    std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
  }
  if (Verbosity() > 3)
  {
    std::cout << "number of duplicates : " << n_dupli << std::endl;
  }
  return iters_out;
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

  static const int _PRINT_THRESHOLD = 2; 

  t_process->restart();


  keyPtrArr ptsPerLayer;
  ptsPerLayer = GetKeyPoints();
  print_reset_time(t_process, Verbosity(), _PRINT_THRESHOLD, "create key-geometry points");
  if (Verbosity()>0) {
    std::cout << " n clusters: " << _cluster_pts.size() << std::endl;
  }

  // get the new links
  linkIterList start_links;
  linkIterArr  links; // array of lists, one for each layer
  std::tie(start_links, links) = CreateBilinks(ptsPerLayer);
  print_reset_time(t_process, Verbosity(), _PRINT_THRESHOLD, "created bilinks");
  if (Verbosity()>0) {
    std::cout << " n start links: " << start_links.size() << std::endl;
    int ntot=0;
    for (const auto& arr : links) { ntot+=arr.size(); }
    std::cout << " n body links: " << ntot << std::endl;
  }

  /*
     Note: The following code segments could very likely use a new algorithm
     Current alorithm:
       Input: links of [layer+1]->[layer]
            1. startlinks: vector<links>: 
               links in which the top layer clusters are not pointed to by 
               any layer above
            2. bodylinks: array<vector<links>>: a sets of links in each layer
               which *do* have an above layer link pointed to them.
               
        Process:
           Start with each starting set of three links (A->B->C) going downward,
           check that passes a curvature test. Use as "starting seed".

           For each starting seed, grow it by links in the following layer.
           Could use a map or make a double-linked list, but don't. We just
           iterate and check.

           In any case, if the seed has multiple matches we just ghost it.
  */

  keyPtrLists seeds  = MakeSeedTripletHeads(start_links, links);
  print_reset_time(t_process, Verbosity(), _PRINT_THRESHOLD, "made triplet heads");
  if (Verbosity()>0) {
    std::cout << " number of seed triplet heads: " << seeds.size() << std::endl;
  }

  GrowSeeds(seeds, links);
  print_reset_time(t_process, Verbosity(), _PRINT_THRESHOLD, "grown seeds");

  auto v2_seeds = RemoveBadClusters(seeds);
  print_reset_time(t_process, Verbosity(), _PRINT_THRESHOLD, "removed bad clusters");
  if (Verbosity()>0) {
    std::cout << " number of seeds: " << v2_seeds.size() << std::endl;
  }

  PublishSeeds(v2_seeds);
  print_reset_time(t_process, Verbosity(), _PRINT_THRESHOLD, "published seeds");

  unsigned int numberofseeds = seeds.size();

  if (Verbosity() > 0)
  {
    std::cout << "number of seeds " << numberofseeds << std::endl;
  }
  if (Verbosity() > 0)
  {
    std::cout << "PHCASeeding Time: " << t_process->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

std::pair<PHCASeeding::linkIterList, PHCASeeding::linkIterArr> PHCASeeding::CreateBilinks(PHCASeeding::keyPtrArr& key_arr)
{

  linkIterList  startLinks; // bilinks at start of chains
  linkIterArr   bodyLinks;  //  bilinks to build chains

  unsigned int sum_size=0;
  for (unsigned int i=0;i<bodyLinks.size();++i) {
    const auto i_size = key_arr[i].size();
    bodyLinks[i].reserve(i_size);
    sum_size += i_size;
  }
  startLinks.reserve(sum_size / bodyLinks.size());

  double cluster_find_time = 0;
  double rtree_query_time = 0;
  double transform_time = 0;
  double compute_best_angle_time = 0;
  double set_insert_time = 0;

  // there are three coord_array (only the current layer is used at a time,
  // but it is filled the same time as the _rtrees_old, which are used two at
  // a time -- the prior padplane row and the next padplain row
  std::array<keyPtrList, 3> iter_arr;
  std::array<std::unordered_set<linkIter>, 2> previous_downlinks_arr;
  std::array<std::set<keyPtr>, 2>  bottom_of_bilink_arr;

  // iterate from outer to inner layers
  const int inner_index = _start_layer - _FIRST_LAYER_TPC + 1;
  const int outer_index = _end_layer - _FIRST_LAYER_TPC - 2;

  // fill the current and prior row coord and ttrees for the first iteration
  int _index_above = (outer_index + 1) % 3;
  int _index_current = (outer_index) % 3;

  // fill the booste trees and remove duplicate clusters
  iter_arr[_index_above]   = FillTree(_rtrees[_index_above], key_arr[outer_index + 1]);
  iter_arr[_index_current] = FillTree(_rtrees[_index_current], key_arr[outer_index]);

  for (int layer_index = outer_index; layer_index >= inner_index; --layer_index)
  {
    // these lines of code will rotates through all three _rtree's in the array,
    // where the old lower becomes the new middle, the old middle the new upper,
    // and the old upper drops out and that _rtree is filled with the new lower
    const unsigned int LAYER = layer_index + _FIRST_LAYER_TPC;
    const double dphi_win = dphi_per_layer[LAYER];
    const double dz_win = dZ_per_layer[LAYER];
    const double dphi_win_p1 = dphi_per_layer[LAYER+1];
    const double dz_win_p1 = dZ_per_layer[LAYER+1];
    int index_above = (layer_index + 1) % 3;
    int index_current = (layer_index) % 3;
    int index_below = (layer_index - 1) % 3;

    iter_arr[index_below] = FillTree(_rtrees[index_below], key_arr[layer_index - 1]);

    // NO DUPLICATES FOUND IN iter_arr

    auto& _rtree_above = _rtrees[index_above];
    const keyPtrList& ptrs = iter_arr[index_current];
    auto& _rtree_below = _rtrees[index_below];

    auto& curr_downlinks = previous_downlinks_arr[layer_index % 2];
    auto& last_downlinks = previous_downlinks_arr[(layer_index + 1) % 2];

    auto& curr_bottom_of_bilink = bottom_of_bilink_arr[layer_index % 2];
    auto& last_bottom_of_bilink = bottom_of_bilink_arr[(layer_index + 1) % 2];

    curr_downlinks.clear();
    curr_bottom_of_bilink.clear();

    // For all the clusters in ptrs, find nearest neighbors in the
    // above and below layers and make links
    // Any link to an above node which matches the same clusters
    // on the previous iteration (to a "below node") becomes a "bilink"
    // Check if this bilink links to a prior bilink or not
    for (const auto& p : ptrs)
    {
      LogDebug(" starting cluster:" << std::endl);
      LogDebug(" z: " << p->z << std::endl);
      LogDebug(" phi: " << p->phi << std::endl);

      /* keyPtrList ClustersAbove; */
      /* keyPtrList ClustersBelow; */
      rtreePairList ClustersAbove;
      rtreePairList ClustersBelow;

      QueryTree(_rtree_below, p, dphi_win, dz_win, ClustersBelow);

      FillTupWinLink(_rtree_below, p);

      QueryTree(_rtree_above, p, dphi_win_p1, dz_win_p1, ClustersAbove);

      t_seed->stop();
      rtree_query_time += t_seed->elapsed();
      t_seed->restart();
      LogDebug(" entries in below layer: " << ClustersBelow.size() << std::endl);
      LogDebug(" entries in above layer: " << ClustersAbove.size() << std::endl);

      std::vector<bool> new_dnlinks(ClustersBelow.size(), true);
      // find the three clusters closest to a straight line
      // (by maximizing the cos of the angle between the (delta_z_,delta_phi) vectors)
      // double minSumLengths = 1e9;
      for (const auto& it_above : ClustersAbove) 
      {
        const auto B = keyPointDiff(*p, *it_above.second);
        bool new_uplink = true;

        for (unsigned int i_below = 0; i_below<ClustersBelow.size(); ++i_below) {
          if (!new_dnlinks[i_below] && !new_uplink) { continue; }

          const auto& it_below = ClustersBelow[i_below].second;
          const auto A = keyPointDiff(*p, *it_below);

          const double A_len_sq = (A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
          const double B_len_sq = (B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
          const double dot_prod = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
          const double cos_angle_sq = dot_prod * dot_prod / A_len_sq / B_len_sq;  
          // also same as cos(angle), where angle is between two vectors

          FillTupWinCosAngle(it_above.second, p, it_below, cos_angle_sq, (dot_prod < 0.));

          constexpr double maxCosPlaneAngle = -0.95;
          constexpr double maxCosPlaneAngle_sq = maxCosPlaneAngle * maxCosPlaneAngle;
          if ((dot_prod < 0.) && (cos_angle_sq > maxCosPlaneAngle_sq))
          {
            if (new_dnlinks[i_below]) {
              new_dnlinks[i_below] = false;
              curr_downlinks.insert({p, it_below}); // std::set, ignores doubles
            }
            if (new_uplink) {
              new_uplink = false;
              linkIter uplink = std::make_pair(it_above.second, p);
              if (last_downlinks.find(uplink) != last_downlinks.end())
              {
                // This is a new bilink
                curr_bottom_of_bilink.insert(p);
                fill_tuple(_tupclus_bilinks, 0, it_above.second);
                fill_tuple(_tupclus_bilinks, 1, p);

                if (last_bottom_of_bilink.find(it_above.second) == last_bottom_of_bilink.end())
                {
                  startLinks.emplace_back(uplink);
                }
                else
                {
                  bodyLinks[layer_index + 1].emplace_back(uplink);
                }
              }
            }  // end new uplink
          } // end check over triplet (cos-angle)
        } // end loop over below matched clusters
      } // end loop over above matched clusters
    } // end loop over all clusters in current layer

    // NOTE:
    // There was some old commented-out code here for allowing layers to be skipped. This
    // may be useful in the future. This chunk of code has been moved towards the
    // end fo the file under the title: "---OLD CODE 0: SKIP_LAYERS---"

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

  // future optimization?
  // Can sort the body links per layer so that links can be binary-searched per layer
  // for now this doesn't seem any faster than just iterating and checking each link,
  // but may become more important in Au+Au data
  /* for (auto& layer : bodyLinks) { std::sort(layer.begin(), layer.end()); } */
  return std::make_pair(startLinks, bodyLinks);
}

double PHCASeeding::getMengerCurvature(PHCASeeding::keyPtr a, PHCASeeding::keyPtr b, PHCASeeding::keyPtr c) const
{
  // Menger curvature = 1/R for circumcircle of triangle formed by most recent three clusters
  // We use here 1/R = 2*sin(breaking angle)/(hypotenuse of triangle)
  double hypot_length = sqrt(square<double>(c->x - a->x) + square<double>(c->y - a->y) + square<double>(c->z - a->z));
  double break_angle = breaking_angle(
      a->x - b->x,
      a->y - b->y,
      a->z - b->z,
      c->x - b->x,
      c->y - b->y,
      c->z - b->z);
  return 2 * sin(break_angle) / hypot_length;
}

PHCASeeding::keyPtrLists PHCASeeding::MakeSeedTripletHeads(
  const PHCASeeding::linkIterList& start_links, 
  const PHCASeeding::linkIterArr& bilinks) const
{
  // form all possible starting 3-cluster tracks (we need that to calculate curvature)
  keyPtrLists seeds;

  // fixme! this is a combinatorial problem -- all copys of triplets will make
  // *many* duplicate tracks from split clusters
  for (auto& startLink : start_links)
  {
    const auto trackHead = startLink.second;
    unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead->key) - _FIRST_LAYER_TPC;
    // the following call with get iterators to all bilinks which match the head
    for (const auto& matchlink : bilinks[trackHead_layer])
    {
      if (matchlink.first == trackHead) {
        seeds.push_back({startLink.first, startLink.second, matchlink.second});
        /* keyPtrList trackSeedTriplet; */
        /* trackSeedTriplet.push_back(startLink.first); */
        /* trackSeedTriplet.push_back(startLink.second); */
        /* trackSeedTriplet.push_back(matchlink.second); */
        /* seeds.push_back(trackSeedTriplet); */

        fill_tuple(_tupclus_seeds, 0, startLink.first);
        fill_tuple(_tupclus_seeds, 1, startLink.second);
        fill_tuple(_tupclus_seeds, 2, matchlink.second);
      }
    }
  }
  return seeds;
}

void PHCASeeding::GrowSeeds(PHCASeeding::keyPtrLists& seeds, const PHCASeeding::linkIterArr& bilinks)
{
  int nsplit_chains = -1;
  unsigned int seed_index = 0;
  while (seed_index < seeds.size()) {
    keyPtrList* seed = &seeds[seed_index];
    seed_index++;

    Checker_dphidz clus_checker(_clusadd_delta_dzdr_window, _clusadd_delta_dphidr2_window, *seed);
    keyPtrList head_keys = {seed->back()};

    while (true) // iterate until break
    {
      // Get all bilinks which fit to the head of the chain
      unsigned int layer = TrkrDefs::getLayer(head_keys[0]->key) - _FIRST_LAYER_TPC;
      keyPtrSet matching_keys{};
      for (const auto& head_key : head_keys)
      {
        // also possible to sort the links and use a sorted search
        // for (auto link = matched_links.first; link != matched_links.second; ++link)
        for (auto& link : bilinks[layer])
        {  // layer for "Index of Layer"
          if (link.first == head_key)
          {
            matching_keys.insert(link.second);
          }
        }
      } // end loop over all bilinks in new layer

      const unsigned int nmatched = matching_keys.size();
      if (nmatched == 1) { // one matched key
        const auto& new_cluster = *matching_keys.begin();
        if (clus_checker.check_cluster(new_cluster)) {
          // first key is a good key
          seed->emplace_back(new_cluster);
          if (seed->size()>=_max_clusters_per_seed) { break; }
          clus_checker.add_cluster();
          head_keys = {new_cluster};
          continue;
        } else { // one matched cluster -- but is not good
          break;
        }
      } else if (nmatched==0) {
        break;
      } 

      // there are multiple matched links
      //   find out which ones are good
      keyPtrList passing_keys{};
      for (const auto& p : matching_keys)
      {  
        // see if the link passes the growth cuts
        if (_split_seeds)
        { FillTupWinGrowSeed(*seed, p); }

        if (clus_checker.check_cluster(p)) 
        { passing_keys.emplace_back(p); } 

        if (_split_seeds)
        { fill_split_chains(*seed, passing_keys, nsplit_chains); }
      } 

      const unsigned int npass = passing_keys.size();
      if (npass > 1) { // most likely case
        if (_split_seeds) { // make new chains up to this link
          for (unsigned int i = 1; i < passing_keys.size(); ++i)
          {
            keyPtrList newseed = {seed->begin(), seed->end()};
            newseed.emplace_back(passing_keys[i]);
            seeds.emplace_back(std::move(newseed));
            seed = &seeds[seed_index]; // seeds might have rearranged
          }
          seed->emplace_back(passing_keys[0]);
          if (seed->size() >= _max_clusters_per_seed) { break; }
          clus_checker.add_cluster(passing_keys[0]);
          head_keys = {passing_keys[0]};
          continue;
        } else {
          // instead of splitting track, add averaged matched cluster position
          clus_checker.add_clusters(passing_keys);
          head_keys = std::move(passing_keys);
          continue;
        }
      } else if (npass == 0) {
        break;
      } else { // there is a single passing seed
        seed->emplace_back(passing_keys[0]);
        if (seed->size() >= _max_clusters_per_seed) { break; }
        clus_checker.add_cluster(passing_keys[0]);
        head_keys = {passing_keys[0]};
        continue;
      }
    } // end growing single seed
         
    if (seed->size() >= _min_clusters_per_seed)
    {
      fill_tuple_with_seed(_tupclus_grown_seeds, *seed);
    }
  }  // end of loop over seeds
}

std::vector<TrackSeed_v2> PHCASeeding::RemoveBadClusters(const PHCASeeding::keyPtrLists& chains) const
{
  if (Verbosity() > 0)
  {
    std::cout << "removing bad clusters" << std::endl;
  }
  std::vector<TrackSeed_v2> clean_chains;

  for (const auto& chain : chains)
  {
    if (chain.size() < 3 || chain.size() < _min_clusters_per_track)
    {
      continue;
    }
    if (Verbosity() > 3)
    {
      std::cout << "chain size: " << chain.size() << std::endl;
    }

    using pvec_t = TrackFitUtils::position_vector_t;

    pvec_t xy_pts;
    std::transform(chain.begin(), chain.end(), std::back_inserter(xy_pts),
      [](const keyPtr& p) -> pvec_t::value_type { return {p->x, p->y}; });

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
    for (const auto& p : chain)
    {
      trackseed.insert_cluster_key(p->key);
    }
    clean_chains.push_back(trackseed);
    if (Verbosity() > 2)
    {
      std::cout << "pushed clean chain with " << trackseed.size_cluster_keys() << " clusters" << std::endl;
    }
  }

  return clean_chains;
}

void PHCASeeding::PublishSeeds(std::vector<TrackSeed_v2>& seeds)
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

  t_process = std::make_unique<PHTimer>("t_seed");
  t_process->stop();

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
  /* _tupclus_links = new TNtuple("links", "links", "event:layer:updown01:x:y:z:delta_z:delta_phi"); */
  _tupclus_bilinks = new TNtuple("bilinks", "bilinks", "event:layer:topbot01:x:y:z");
  _tupclus_seeds = new TNtuple("seeds", "3 bilink seeds cores", "event:layer:seed012:x:y:z");
  _tupclus_grown_seeds = new TNtuple("grown_seeds", "grown seeds", "event:layer:seednum05:x:y:z:dzdr:d2phidr2");
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

void PHCASeeding::fill_split_chains(const PHCASeeding::keyPtrList& seed, const PHCASeeding::keyPtrList& add_links, int& n_tupchains) const
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
      r2 = r1;
      phi2 = phi1;
    }
    if (has_0)
    {
      has_1 = true;
      z1 = z0;
      r1 = r0;
      phi1 = phi0;
    }
    /* auto link_pos = pos.at(link); */
    has_0 = true;
    x0 = link->x;
    y0 = link->y;
    z0 = link->z;
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
    float dphi01 = wrap_dphi(phi0, phi1);
    float dphi12 = wrap_dphi(phi1, phi2);
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
    /* auto link_pos = pos.at(link); */
    float xt = link->x;
    float yt = link->y;
    float zt = link->z;
    float phit = atan2(yt, xt);
    float rt = sqrt(xt * xt + yt * yt);
    float dr_t0 = rt - r0;

    float dzdr = (zt - z0) / (rt - r0);
    float dphit0 = wrap_dphi(phi0, phit);

    float d2phidr2 = dphit0 / dr_t0 / dr_t0 - dphidr01;
    _tup_chainfork->Fill(_tupout_count, n_tupchains, TrkrDefs::getLayer(link), xt, yt, zt, dzdr, d2phidr2, index, nlinks);
  }
}

#else
void PHCASeeding::fill_split_chains(const PHCASeeding::keyPtrList& /*chain*/, const PHCASeeding::keyPtrList& /*links*/, int& /*nchains*/) const {};
#endif

#if defined(_PHCASEEDING_CLUSTERLOG_TUPOUT_)
void PHCASeeding::write_tuples()
{
  _f_clustering_process->cd();
  _tupclus_all->Write();
  /* _tupclus_links->Write(); */
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

void PHCASeeding::fill_tuple(TNtuple* tup, float val, PHCASeeding::keyPtr p) const
{
  tup->Fill(_tupout_count, TrkrDefs::getLayer(p->key), val, p->x, p->y, p->z);
}

void PHCASeeding::fill_tuple_with_seed(TNtuple* tup, const PHCASeeding::keyPtrList& seed) const
{
  keyPtrList trio {seed.begin(), seed.begin()+3};
  Checker_dphidz checker(_clusadd_delta_dzdr_window, _clusadd_delta_dphidr2_window, trio);

  int cnt = 0;
  // also output the dR/dZ and dphi2/dR2
  for (unsigned int i=0; i<seed.size(); ++i) {
    if (i>2) {
      checker.add_cluster(seed[i]);
    }
    tup->Fill(_tupout_count, TrkrDefs::getLayer(seed[i]->key), cnt, seed[i]->x, seed[i]->y, seed[i]->z, checker.calc_dzdr(), checker.calc_d2phidr2());
    ++cnt;
  }
}

void PHCASeeding::process_tupout_count()
{
  _tupout_count += 1;
  if (_tupout_count != 0) return;
  _search_windows->Fill(_neighbor_z_width, _neighbor_phi_width, _start_layer, _end_layer, _clusadd_delta_dzdr_window, _clusadd_delta_dphidr2_window);
}

void PHCASeeding::FillTupWinGrowSeed(const PHCASeeding::keyPtrList& seed, const PHCASeeding::keyPtr& p) const
{
  keyPtr head = seed.back();
  keyPtr prev = seed.rbegin()[1];
  float x1 = head->x;
  float y1 = head->y;
  float z1 = head->z;
  float x2 = prev->x;
  float y2 = prev->y;
  float z2 = prev->z;
  float dr_12 = sqrt(x1 * x1 + y1 * y1) - sqrt(x2 * x2 + y2 * y2);
  /* TrkrDefs::cluskey testCluster = link.second; */
  float xt = p->x;
  float yt = p->y;
  float zt = p->z;
  float dr_t1 = sqrt(xt * xt + yt * yt) - sqrt(x1 * x1 + y1 * y1);
  float dzdr_12 = (z1 - z2) / dr_12;
  float dzdr_t1 = (zt - z1) / dr_t1;
  // if (fabs(dzdr_12 - dzdr_t1) > _clusadd_delta_dzdr_window)) // then fail this link

  keyPtr third_pos = seed.rbegin()[2];
  float x3 = third_pos->x;
  float y3 = third_pos->y;
  float z3 = third_pos->z;
  float dr_23 = sqrt(x2 * x2 + y2 * y2) - sqrt(x3 * x3 + y3 * y3);
  float phi1 = atan2(y1, x1);
  float phi2 = atan2(y2, x2);
  float phi3 = atan2(y3, x3);
  float dphi12 = wrap_dphi(phi2, phi1);
  float dphi23 = wrap_dphi(phi3, phi2); 
  float d2phidr2_123 = dphi12 / (dr_12 * dr_12) - dphi23 / (dr_23 * dr_23);
  float dphit1 = std::fmod(atan2(yt, xt) - atan2(y1, x1), M_PI);
  float d2phidr2_t12 = dphit1 / (dr_t1 * dr_t1) - dphi12 / (dr_12 * dr_12);
  _tupwin_seed23->Fill(_tupout_count,
                       (TrkrDefs::getLayer(seed.rbegin()[1]->key)), x2, y2, z2,
                       (TrkrDefs::getLayer(seed.rbegin()[2]->key)), x3, y3, z3);
  _tupwin_seedL1->Fill(_tupout_count,
                       (TrkrDefs::getLayer(p->key)), xt, yt, zt, // ? not sure what orig
                       (TrkrDefs::getLayer(seed.back()->key)), x1, y1, z1,
                       dzdr_12, dzdr_t1, fabs(dzdr_12 - dzdr_t1),
                       d2phidr2_123, d2phidr2_t12, fabs(d2phidr2_123 - d2phidr2_t12));
}

void PHCASeeding::FillTupWinLink(PHCASeeding::boost_rtree& rtree_below, const PHCASeeding::keyPtr p) const
{
  double StartPhi = p->phi;
  double StartZ = p->z;
  // Fill TNTuple _tupwin_link
  rtreePairList ClustersBelow;
  QueryTree(rtree_below, p, 1., 20, ClustersBelow);

  for (const auto& pkey : ClustersBelow)
  {
    const auto new_phi = pkey.second->phi;
    double dphi = wrap_dphi(StartPhi, new_phi);
    double dZ = pkey.second->z - StartZ;
    _tupwin_link->Fill(_tupout_count, TrkrDefs::getLayer(p->key), p->x,p->y,p->z, TrkrDefs::getLayer(pkey.second->key), pkey.second->x,pkey.second->y,pkey.second->z, dphi, dZ);
  }
}

void PHCASeeding::FillTupWinCosAngle(const PHCASeeding::keyPtr A, const PHCASeeding::keyPtr B, const PHCASeeding::keyPtr C, double cos_angle_sq, bool isneg) const
{
  // A is top cluster, B the middle, C the bottom
  // a,b,c are the positions

  _tupwin_cos_angle->Fill(_tupout_count,
                          TrkrDefs::getLayer(A->key), A->x, A->y, A->z,
                          TrkrDefs::getLayer(B->key), B->x, B->y, B->z,
                          TrkrDefs::getLayer(C->key), C->x, C->y, C->z,
                          (isneg ? -1 : 1) * sqrt(cos_angle_sq));
}
#else
void PHCASeeding::write_tuples(){};
void PHCASeeding::fill_tuple(TNtuple* /**/, float /**/, PHCASeeding::keyPtr /**/) const {};
void PHCASeeding::fill_tuple_with_seed(TNtuple* /**/, const PHCASeeding::keyPtrList& /**/) const {};
void PHCASeeding::process_tupout_count(){};
void PHCASeeding::FillTupWinGrowSeed(const PHCASeeding::keyPtrList& /**/, const PHCASeeding::keyPtr& /**/) const {};
void PHCASeeding::FillTupWinLink(PHCASeeding::boost_rtree& /**/, const PHCASeeding::keyPtr /**/) const {};
void PHCASeeding::FillTupWinCosAngle(PHCASeeding::keyPtr /**/, PHCASeeding::keyPtr /**/, PHCASeeding::keyPtr /**/, double /**/, bool /**/) const {}; 
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
    QueryTree_old(_rtree,
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
      QueryTree_old(_rtree,
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
