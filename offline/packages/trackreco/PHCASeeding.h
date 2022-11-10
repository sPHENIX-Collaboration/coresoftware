#ifndef TRACKRECO_PHCASEEDING_H
#define TRACKRECO_PHCASEEDING_H

/*!
 *  \file PHCASeeding.cc
 *  \brief Track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail 
 *  \author Michael Peters & Christof Roland
 */

//begin

#include "PHTrackSeeding.h"      // for PHTrackSeeding
#include "ALICEKF.h"

#include <tpc/TpcDistortionCorrection.h>

#include <trackbase/TrkrDefs.h>  // for cluskey
#include <trackbase/ActsGeometry.h>

#include <phool/PHTimer.h>  // for PHTimer

#include <Eigen/Core>
#include <Eigen/Dense>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/geometry/geometries/box.hpp>    // for box
#pragma GCC diagnostic pop

#include <boost/geometry/geometries/point.hpp>  // for point
#include <boost/geometry/index/rtree.hpp>       // for ca

#include <cmath>     // for M_PI
#include <cstdint>  // for uint64_t
#include <map>       // for map
#include <memory>
#include <set>
#include <string>    // for string
#include <utility>   // for pair
#include <unordered_set>
#include <vector>    // for vector

class PHCompositeNode;  
class PHTimer;
class SvtxTrack_v3;
class TpcDistortionCorrectionContainer;
class TrkrCluster;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using point = bg::model::point<float, 3, bg::cs::cartesian>;
using box = bg::model::box<point>;
using pointKey = std::pair<point, TrkrDefs::cluskey>;
using coordKey = std::pair<std::array<float,3>, TrkrDefs::cluskey>;
using keylink = std::array<coordKey,2>;
using keylist = std::vector<TrkrDefs::cluskey>;
using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;

class PHCASeeding : public PHTrackSeeding
{
 public:
  PHCASeeding(
      const std::string &name = "PHCASeeding",
      unsigned int start_layer = 7,
      unsigned int end_layer = 55,
      unsigned int min_nhits_per_cluster = 0,
      unsigned int min_clusters_per_track = 5,
      const unsigned int nlayers_maps = 3,
      const unsigned int nlayers_intt = 4,
      const unsigned int nlayers_tpc = 48,
      float neighbor_phi_width = .02,
      float neighbor_eta_width = .01,
      float maxSinPhi = 0.999,
      float cosTheta_limit = -0.8);

  ~PHCASeeding() override {}
  void SetLayerRange(unsigned int layer_low, unsigned int layer_up) {_start_layer = layer_low; _end_layer = layer_up;}
  void SetSearchWindow(float eta_width, float phi_width) {_neighbor_eta_width = eta_width; _neighbor_phi_width = phi_width;}
  void SetMinHitsPerCluster(unsigned int minHits) {_min_nhits_per_cluster = minHits;}
  void SetMinClustersPerTrack(unsigned int minClus) {_min_clusters_per_track = minClus;}

  void set_field_dir(const double rescale)
  {
    std::cout << "rescale: " << rescale << std::endl;
    _fieldDir = 1;
    if(rescale > 0)
      _fieldDir = -1;     
  }

  void useConstBField(bool opt){_use_const_field = opt;}
  void useFixedClusterError(bool opt){_use_fixed_clus_err = opt;}
  void setFixedClusterError(int i, double val){_fixed_clus_err.at(i) = val;}

 protected:
  int Setup(PHCompositeNode *topNode) override;
  int Process(PHCompositeNode *topNode) override;
  int InitializeGeometry(PHCompositeNode *topNode);
  int FindSeedsLayerSkip(double cosTheta_limit);
  int End() override;

 private:
  
  enum skip_layers {on, off};
  
  /// tpc distortion correction utility class
  TpcDistortionCorrection m_distortionCorrection;

  /// get global position for a given cluster
  /**
   * uses ActsTransformation to convert cluster local position into global coordinates
   * incorporates TPC distortion correction, if present
   */
  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey, TrkrCluster*) const;

  PositionMap FillTree();
  int FindSeedsWithMerger(const PositionMap&);
  std::pair<std::vector<std::unordered_set<keylink>>,std::vector<std::unordered_set<keylink>>> CreateLinks(const std::vector<coordKey>& clusters, const PositionMap& globalPositions) const;
  std::vector<std::vector<keylink>> FindBiLinks(const std::vector<std::unordered_set<keylink>>& belowLinks, const std::vector<std::unordered_set<keylink>>& aboveLinks) const;
  std::vector<keylist> FollowBiLinks(const std::vector<std::vector<keylink>>& bidirectionalLinks, const PositionMap& globalPositions) const;
  void QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax, std::vector<pointKey> &returned_values) const;
  std::vector<TrackSeed_v1> RemoveBadClusters(const std::vector<keylist>& seeds, const PositionMap& globalPositions) const;
  double getMengerCurvature(TrkrDefs::cluskey a, TrkrDefs::cluskey b, TrkrDefs::cluskey c, const PositionMap& globalPositions) const;
  
  void publishSeeds(const std::vector<TrackSeed_v1>& seeds);

  //int _nlayers_all;
  //unsigned int _nlayers_seeding;
  //std::vector<int> _seeding_layer;

  const unsigned int _nlayers_maps;
  const unsigned int _nlayers_intt;
  const unsigned int _nlayers_tpc;
  unsigned int _start_layer;
  unsigned int _end_layer;
  unsigned int _min_nhits_per_cluster;
  unsigned int _min_clusters_per_track;
//  float _cluster_z_error;
//  float _cluster_alice_y_error;
  float _neighbor_phi_width;
  float _neighbor_eta_width;
  float _max_sin_phi;
  float _cosTheta_limit;
  double _rz_outlier_threshold = 0.1;
  double _xy_outlier_threshold = 0.1;
  double _fieldDir = -1;
  bool _use_const_field = false;
  bool _use_fixed_clus_err = false;
  std::array<double,3> _fixed_clus_err = {.1,.1,.1};

  /// acts geometry
  ActsGeometry *tGeometry{nullptr};

  /// distortion correction container
  TpcDistortionCorrectionContainer* m_dcc = nullptr;
  
  std::unique_ptr<ALICEKF> fitter;
 
  std::unique_ptr<PHTimer> t_seed;
  std::unique_ptr<PHTimer> t_fill;
  bgi::rtree<pointKey, bgi::quadratic<16>> _rtree;
};

#endif
