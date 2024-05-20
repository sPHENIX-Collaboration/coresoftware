#ifndef TRACKRECO_PHCASEEDING_H
#define TRACKRECO_PHCASEEDING_H

/*!
 *  \file PHCASeeding.cc
 *  \brief Track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail
 *  \author Michael Peters & Christof Roland
 */

// Statements for if we want to save out the intermediary clustering steps
/* #define _CLUSTER_LOG_TUPOUT_ */
#define _PHCASEEDING_TIMER_OUT_

#include "ALICEKF.h"
#include "PHTrackSeeding.h"  // for PHTrackSeeding

#include <tpc/TpcDistortionCorrection.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>  // for cluskey

#include <phool/PHTimer.h>  // for PHTimer

#include <Eigen/Core>
#include <Eigen/Dense>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/geometry/geometries/box.hpp>  // for box
#pragma GCC diagnostic pop

#include <boost/geometry/geometries/point.hpp>  // for point
#include <boost/geometry/index/rtree.hpp>       // for ca

#include <cmath>    // for M_PI
#include <cstdint>  // for uint64_t
#include <map>      // for map
#include <unordered_map>      // for map
#include <memory>
#include <set>
#include <string>  // for string
#include <unordered_set>
#include <utility>  // for pair
#include <vector>   // for vector

#if defined(_CLUSTER_LOG_TUPOUT_)
#include <TFile.h>
#include <TNtuple.h>
#endif

class PHCompositeNode;
class PHTimer;
class SvtxTrack_v3;
class TpcDistortionCorrectionContainer;
class TrkrCluster;


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


class PHCASeeding : public PHTrackSeeding
{
 public:
   static const int _NLAYERS_TPC = 48;
   static const int _FIRST_LAYER_TPC = 7;
   // move `using` statements inside of the class to avoid polluting the global namespace

   using point = bg::model::point<float, 2, bg::cs::cartesian>;
   using box = bg::model::box<point>;
   using pointKey = std::pair<point, TrkrDefs::cluskey>;
   using coordKey = std::pair<std::array<float, 2>, TrkrDefs::cluskey>; // just use phi and eta, no longer needs the layer

   using keyList = std::vector<TrkrDefs::cluskey>;
   using keyLists = std::vector<keyList>;
   using keyListPerLayer = std::array<keyList, _NLAYERS_TPC>;

   using keyLink = std::pair<TrkrDefs::cluskey, TrkrDefs::cluskey>;
   using keyLinks = std::vector<keyLink>;
   using keyLinkPerLayer = std::array<std::vector<keyLink>,_NLAYERS_TPC>;

   using PositionMap = std::unordered_map<TrkrDefs::cluskey, Acts::Vector3>;

   PHCASeeding(
      const std::string& name = "PHCASeeding",
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
  void SetLayerRange(unsigned int layer_low, unsigned int layer_up)
  {
    _start_layer = layer_low;
    _end_layer = layer_up;
  }
  void SetSearchWindow(float eta_width, float phi_width)
  {
    _neighbor_eta_width = eta_width;
    _neighbor_phi_width = phi_width;
  }
  void SetClusAdd_delta_window(float _dzdr_cutoff, float _dphidr2_cutoff) 
  {
    _clusadd_delta_dzdr_window = _dzdr_cutoff;
    _clusadd_delta_dphidr2_window = _dphidr2_cutoff;
  }
  void SetMinHitsPerCluster(unsigned int minHits) { _min_nhits_per_cluster = minHits; }
  void SetMinClustersPerTrack(unsigned int minClus) { _min_clusters_per_track = minClus; }
  void SetNClustersPerSeedRange(unsigned int minClus, unsigned int maxClus) 
  { _min_clusters_per_seed = minClus;  _max_clusters_per_seed = maxClus; 
    if (_min_clusters_per_seed < 3) {
      std::cout << " Error in SetNClustersPerSeedRange: " << __FILE__
        << " min value cannot be less than three." << std::endl;
      assert(false);
    }
  }

  void set_field_dir(const double rescale)
  {
    std::cout << "PHCASeeding::set_field_dir rescale: " << rescale << std::endl;
    _fieldDir = 1;
    if (rescale > 0)
      _fieldDir = -1;
  }

  void magFieldFile(const std::string& fname) { m_magField = fname; }
  void useConstBField(bool opt) { _use_const_field = opt; }
  void constBField(float b) { _const_field = b; }
  void useFixedClusterError(bool opt) { _use_fixed_clus_err = opt; }
  void setFixedClusterError(int i, double val) { _fixed_clus_err.at(i) = val; }
  void set_pp_mode(bool mode) { _pp_mode = mode; }

 protected:
  int Setup(PHCompositeNode* topNode) override;
  int Process(PHCompositeNode* topNode) override;
  int InitializeGeometry(PHCompositeNode* topNode);
  int FindSeedsLayerSkip(double cosTheta_limit);
  int End() override;

 private:
  bool _save_clus_proc = false;

#if defined(_CLUSTER_LOG_TUPOUT_)
  TFile* _f_clustering_process = nullptr;
  int      _tupout_count=-1;
  // Save the steps of the clustering
  TNtuple* _tupclus_all = nullptr;    // all clusters
  TNtuple* _tupclus_links = nullptr; //  clusters which are linked
  TNtuple* _tupclus_bilinks = nullptr; // bi-linked clusters
  TNtuple* _tupclus_seeds = nullptr; // seed (outermost three bi-linked chain
  TNtuple* _tupclus_grown_seeds = nullptr; // seeds saved out
#endif

  // have a comparator to search vector of sorted bilinks for links with starting
  // a key of a given value
  class CompKeyToBilink {
   public:
    bool operator()(const keyLink& p, const TrkrDefs::cluskey& val) const
    {
      return p.first < val;
    }
    bool operator()(const TrkrDefs::cluskey& val, const keyLink& p) const
    {
      return val < p.first;
    }
  };
  std::pair<std::vector<keyLink>::iterator, std::vector<keyLink>::iterator> FindBilinks(const TrkrDefs::cluskey& key);

  /// tpc distortion correction utility class
  TpcDistortionCorrection m_distortionCorrection;

  /// get global position for a given cluster
  /**
   * uses ActsTransformation to convert cluster local position into global coordinates
   * incorporates TPC distortion correction, if present
   */
  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey, TrkrCluster*) const;
  std::pair<PositionMap, keyListPerLayer> FillGlobalPositions();
  std::pair<keyLinks, keyLinkPerLayer> CreateBiLinks(const PositionMap& globalPositions, const keyListPerLayer& ckeys);
  PHCASeeding::keyLists FollowBiLinks( const keyLinks& trackSeedPairs, const keyLinkPerLayer& bilinks, const PositionMap& globalPositions) const;
  std::vector<coordKey> FillTree(bgi::rtree<pointKey,bgi::quadratic<16>>&, const keyList&, const PositionMap&, int layer);
  int FindSeedsWithMerger(const PositionMap&, const keyListPerLayer&);

  void QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>>& rtree, double phimin, double etamin, double phimax, double etamax, std::vector<pointKey>& returned_values) const;
  std::vector<TrackSeed_v2> RemoveBadClusters(const std::vector<keyList>& seeds, const PositionMap& globalPositions) const;
  double getMengerCurvature(TrkrDefs::cluskey a, TrkrDefs::cluskey b, TrkrDefs::cluskey c, const PositionMap& globalPositions) const;

  void publishSeeds(const std::vector<TrackSeed_v2>& seeds);

  // int _nlayers_all;
  // unsigned int _nlayers_seeding;
  // std::vector<int> _seeding_layer;

  const unsigned int _nlayers_maps;
  const unsigned int _nlayers_intt;
  const unsigned int _nlayers_tpc;
  unsigned int _start_layer;
  unsigned int _end_layer;
  unsigned int _min_nhits_per_cluster;
  unsigned int _min_clusters_per_track;
  unsigned int _max_clusters_per_seed = 6;
  unsigned int _min_clusters_per_seed = 6;
  //  float _cluster_z_error;
  //  float _cluster_alice_y_error;
  float _neighbor_phi_width;
  float _neighbor_eta_width;
  float _clusadd_delta_dzdr_window = 0.5;
  float _clusadd_delta_dphidr2_window = 0.005;
  float _max_sin_phi;
  float _cosTheta_limit;
  double _rz_outlier_threshold = 0.1;
  double _xy_outlier_threshold = 0.1;
  double _fieldDir = -1;
  bool _use_const_field = false;
  float _const_field = 1.4;
  bool _use_fixed_clus_err = false;
  bool _pp_mode = false;
  std::array<double, 3> _fixed_clus_err = {.1, .1, .1};

  std::string m_magField;

  /// acts geometry
  ActsGeometry* tGeometry{nullptr};

  /// distortion correction container
  TpcDistortionCorrectionContainer* m_dcc = nullptr;

  std::unique_ptr<ALICEKF> fitter;

  std::unique_ptr<PHTimer> t_seed;
  std::unique_ptr<PHTimer> t_fill;
  std::unique_ptr<PHTimer> t_makebilinks;
  std::unique_ptr<PHTimer> t_makeseeds;
  /* std::array<bgi::rtree<pointKey, bgi::quadratic<16>>, _NLAYERS_TPC> _rtrees; */
  std::array<bgi::rtree<pointKey,bgi::quadratic<16>>, 3> _rtrees; // need three layers at a time
};

#endif
