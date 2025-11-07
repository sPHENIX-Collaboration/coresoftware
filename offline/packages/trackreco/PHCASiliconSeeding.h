#ifndef TRACKRECO_PHCASILICONSEEDING_H
#define TRACKRECO_PHCASILICONSEEDING_H

/*!
 *  \file PHCASiliconSeeding.cc
 *  \brief Silicon track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail
 *  \author Michael Peters
 */

#include "PHTrackSeeding.h"  // for PHTrackSeeding

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterCrossingAssoc.h>
#include <trackbase/TrkrDefs.h>  // for cluskey

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4mvtx/PHG4MvtxMisalignment.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/geometry/geometries/box.hpp>  // for box
#pragma GCC diagnostic pop

#include <boost/geometry/geometries/point.hpp>  // for point
#include <boost/geometry/index/rtree.hpp>       // for ca

#include <cmath>    // for M_PI
#include <cstdint>  // for uint64_t
#include <map>      // for map
#include <memory>
#include <set>
#include <string>  // for string
#include <unordered_map>
#include <utility>  // for pair
#include <vector>   // for vector

class PHCompositeNode;
class TrkrCluster;
class TrackSeed;
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
      float eta_allowance = 1.1);

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
  void SetSearchWindow(float eta_allowance, float phi_width)
  {
    _eta_allowance = eta_allowance;
    _drdz_allowance = (eta_allowance < 1E-5) ? (2. * std::exp(-5.0)) / (1 - std::exp(-1. * 5.0 * 2)) : (2. * std::exp(-eta_allowance)) / (1 - std::exp(-1. * eta_allowance * 2));  // If eta_allowance is very small, set to the value for eta = 5.0 for large allowance
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
  void SetMinMVTXClusters(unsigned int minMVTX)
  {
    _min_mvtx_clusters = minMVTX;
  }
  void SetMinINTTClusters(unsigned int minINTT)
  {
    _min_intt_clusters = minINTT;
  }
  void SetTrackMapName(const std::string& trackmap_name)
  {
    _module_trackmap_name = trackmap_name;
  }
  void Identify() const
  {
    std::cout << "----- Configuration parameters -----" << std::endl;
    std::cout << " - Start layer: " << _start_layer << std::endl;
    std::cout << " - End layer: " << _end_layer << std::endl;
    std::cout << " - MVTX strobe ID range: [" << _lowest_allowed_strobeid << ", " << _highest_allowed_strobeid << "]" << std::endl;
    std::cout << " - Search window eta allowance: " << _eta_allowance << ", drdz allowance: " << _drdz_allowance << std::endl;
    std::cout << " - Search window neighbor phi width: " << _neighbor_phi_width << std::endl;
    std::cout << " - Propagate max DCAxy: " << _propagate_max_dcaxy << std::endl;
    std::cout << " - Propagate max DCAz: " << _propagate_max_dcaz << std::endl;
    std::cout << " - Max cos(triplet breaking angle): " << _max_cos_angle << std::endl;
    std::cout << " - Use best triplet: " << (_use_best ? "true" : "false") << std::endl;
    std::cout << " - Require INTT consistency: " << (_require_INTT_consistency ? "true" : "false") << std::endl;
    std::cout << " - Minimum clusters per track: " << _min_clusters_per_track << std::endl;
    std::cout << " - Minimum MVTX clusters: " << _min_mvtx_clusters << std::endl;
    std::cout << " - Minimum INTT clusters: " << _min_intt_clusters << std::endl;
    std::cout << " - Track Map Name: " << _module_trackmap_name << std::endl;
    std::cout << "------------------------------------" << std::endl;
  }

 protected:
  int Setup(PHCompositeNode* topNode) override;
  int Process(PHCompositeNode* topNode) override;
  int InitializeGeometry(PHCompositeNode* topNode);
  int End() override;

 private:
  unsigned int _start_layer;
  unsigned int _end_layer;
  float _strobe_width = 9.9; // mus
  unsigned int _min_clusters_per_track;
  unsigned int _min_mvtx_clusters = 2;
  unsigned int _min_intt_clusters = 1;

  float _eta_allowance = 1.1;
  float _drdz_allowance = (2. * std::exp(-_eta_allowance)) / (1 - std::exp(-1. * _eta_allowance * 2));  // default for eta allowance of 1.1
  float _neighbor_phi_width;

  int _lowest_allowed_strobeid = -5;
  int _highest_allowed_strobeid = 5;

  float _propagate_max_dcaxy = 0.4;
  float _propagate_max_dcaz = 1.;

  float _max_cos_angle = -0.95;
  bool _use_best = true;
  bool _require_INTT_consistency = true;

  std::array<float, 55> dphi_per_layer{};
  std::array<float, 55> max_dcaxy_perlayer{};
  std::array<float, 55> max_dcaz_perlayer{};
  std::array<float, 7> radius_per_layer{};  // radius of each layer

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
  bool timingMismatch(const TrackSeed& seed) const;
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

  // set up layer radii
  void SetupDefaultLayerRadius();
  float GetMvtxRadiusByPhi(float clusphi, int layer) const;

  /// acts geometry
  ActsGeometry* m_tGeometry{nullptr};
  // cylinder geometry for mvtx and intt
  PHG4CylinderGeomContainer* geom_container_mvtx = nullptr;
  PHG4CylinderGeomContainer* geom_container_intt = nullptr;
  PHG4MvtxMisalignment* mvtxmisalignment = nullptr;
  std::vector<double> v_globaldisplacement = {0., 0., 0.};
  float radius_displacement = 0.;

  //  TrackSeedContainer *m_seedContainer = nullptr;
  TrkrClusterContainer* m_clusterMap = nullptr;
  TrkrClusterCrossingAssoc* m_clusterCrossingMap = nullptr;

  std::string _module_trackmap_name = "SiliconTrackSeedContainer";

  std::vector<bgi::rtree<pointKey, bgi::quadratic<16>>> _rtrees;  // need three layers at a time
};

#endif
