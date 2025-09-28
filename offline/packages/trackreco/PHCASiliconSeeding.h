#ifndef TRACKRECO_PHCASILICONSEEDING_H
#define TRACKRECO_PHCASILICONSEEDING_H

/*!
 *  \file PHCASiliconSeeding.cc
 *  \brief Silicon track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail
 *  \author Michael Peters
 */

#include "PHTrackSeeding.h"  // for PHTrackSeeding

#include <trackbase/TrkrDefs.h>  // for cluskey
#include <trackbase/TrkrClusterCrossingAssoc.h>
#include <trackbase/ActsGeometry.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/geometry/geometries/box.hpp>  // for box
#pragma GCC diagnostic pop

#include <boost/geometry/geometries/point.hpp>  // for point
#include <boost/geometry/index/rtree.hpp>       // for ca

#include <cmath>    // for M_PI
#include <cstdint>  // for uint64_t
#include <map>      // for map
#include <unordered_map>
#include <memory>
#include <set>
#include <string>         // for string
#include <utility>  // for pair
#include <vector>   // for vector

class PHCompositeNode;
class TrkrCluster;
class TrackSeed_v2;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class PHCASiliconSeeding : public PHTrackSeeding
{
 public:
  using point = bg::model::point<float, 2, bg::cs::cartesian>;
  using box = bg::model::box<point>;
  using pointKey = std::pair<point, TrkrDefs::cluskey>;                 // phi and z in the key
  using coordKey = std::pair<std::array<float, 2>, TrkrDefs::cluskey>;  // just use phi and Z, no longer needs the layer

  using keyList = std::vector<TrkrDefs::cluskey>;
  using keyLists = std::vector<keyList>;
  using keyListPerLayer = std::array<keyList, 55>;

  using PositionMap = std::unordered_map<TrkrDefs::cluskey, Acts::Vector3>;

  PHCASiliconSeeding(
      const std::string& name = "PHCASiliconSeeding",
      unsigned int start_layer = 7,
      unsigned int end_layer = 55,
      unsigned int min_clusters_per_track = 5,
      float neighbor_phi_width = .02,
      float neighbor_z_width = .01
  );

  ~PHCASiliconSeeding() override {}
  void SetLayerRange(unsigned int layer_low, unsigned int layer_up)
  {
    _start_layer = layer_low;
    _end_layer = layer_up;
  }
  void SetMVTXStrobeIDRange(int low_strobe, int high_strobe)
  {
    _lowest_allowed_strobeid = low_strobe;
    _highest_allowed_strobeid = high_strobe;
  }
  void SetSearchWindow(float z_width, float phi_width)
  {
    _neighbor_z_width = z_width;
    _neighbor_phi_width = phi_width;
  }
  void SetPropagateMaxDCAxy(float dcaxy)
  {
    _propagate_max_dcaxy = dcaxy;
  }
  void SetPropagateMaxDCAz(float dcaz)
  {
    _propagate_max_dcaz = dcaz;
  }
  void SetMaxCosTripletBreakingAngle(float cos_angle)
  {
    _max_cos_angle = cos_angle;
  }
  void SetAlgoUseBestTriplet(bool use_best)
  {
    _use_best = use_best;
  }
  void SetRequireINTTConsistency(bool req)
  {
    _require_INTT_consistency = req;
  }
  void SetMinClustersPerTrack(unsigned int minClus) 
  { 
    _min_clusters_per_track = minClus; 
  }

  void set_track_map_name(const std::string &map_name) { trackmapname = map_name; }

 protected:
  int Setup(PHCompositeNode* topNode) override;
  int Process(PHCompositeNode* topNode) override;
  int InitializeGeometry(PHCompositeNode* topNode);
  int End() override;

 private:
  unsigned int _start_layer;
  unsigned int _end_layer;
  unsigned int _min_clusters_per_track;
  float _neighbor_phi_width;
  float _neighbor_z_width;

  int _lowest_allowed_strobeid = -5;
  int _highest_allowed_strobeid = 5;

  float _propagate_max_dcaxy = 0.4;
  float _propagate_max_dcaz = 1.;

  float _max_cos_angle = -0.95;
  bool _use_best = true;
  bool _require_INTT_consistency = true;

  std::array<float, 55> dZ_per_layer{};
  std::array<float, 55> dphi_per_layer{};
  std::array<float, 55> max_dcaxy_perlayer{};
  std::array<float, 55> max_dcaz_perlayer{};

  struct Triplet
  {
    TrkrDefs::cluskey bottom;
    TrkrDefs::cluskey center;
    TrkrDefs::cluskey top;
  };

  /// get global position for a given cluster
  /**
   * uses ActsTransformation to convert cluster local position into global coordinates
   * incorporates TPC distortion correction, if present
   */
  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey, TrkrCluster*) const;
  std::pair<PositionMap, keyListPerLayer> FillGlobalPositions();

  std::vector<std::vector<Triplet>> CreateLinks(const PHCASiliconSeeding::PositionMap& globalPositions, const PHCASiliconSeeding::keyListPerLayer& ckeys);
  std::vector<keyList> FollowLinks(const std::vector<std::vector<Triplet>>& triplets);
  std::vector<coordKey> FillTree(bgi::rtree<pointKey, bgi::quadratic<16>>&, const keyList&, const PositionMap&, int layer);
  int FindSeeds(const PositionMap&, const keyListPerLayer&);

  void QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>>& rtree, double phimin, double zmin, double phimax, double zmax, std::vector<pointKey>& returned_values) const;
  float getSeedQuality(const TrackSeed_v2& seed, const PositionMap& globalPositions) const;
  void HelixPropagate(std::vector<TrackSeed_v2>& seeds, const PositionMap& globalPositions) const;
  void FitSeed(TrackSeed_v2& seed, const PositionMap& globalPositions) const;
  std::vector<TrackSeed_v2> FitSeeds(const std::vector<keyList>& seeds, const PositionMap& globalPositions) const;

  std::set<short> GetINTTClusterCrossings(const TrkrDefs::cluskey ckey) const;
  short GetCleanINTTClusterCrossing(const TrkrDefs::cluskey ckey) const;
  int GetClusterTimeIndex(const TrkrDefs::cluskey ckey) const;
  bool ClusterTimesAreCompatible(const uint8_t trkr_id, const int time_index, const TrkrDefs::cluskey ckey) const;
  bool ClusterTimesAreCompatible(const TrkrDefs::cluskey clus_a, const TrkrDefs::cluskey clus_b) const;

  void publishSeeds(const std::vector<TrackSeed_v2>& seeds) const;

  /// acts geometry
  ActsGeometry* m_tGeometry{nullptr};

  std::string trackmapname = "SiliconTrackSeedContainer"; 
  //  TrackSeedContainer *m_seedContainer = nullptr;
  TrkrClusterContainer *m_clusterMap = nullptr;
  TrkrClusterCrossingAssoc *m_clusterCrossingMap = nullptr;

  std::vector<bgi::rtree<pointKey, bgi::quadratic<16>>> _rtrees;  // need three layers at a time
};

#endif
