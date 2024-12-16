#ifndef TRACKRECO_PHCASEEDING_H
#define TRACKRECO_PHCASEEDING_H

/*!
 *  \file PHCASeeding.cc
 *  \brief Track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail
 *  \author Michael Peters & Christof Roland
 */

// Statements for if we want to save out the intermediary clustering steps
/* #define _PHCASEEDING_CLUSTERLOG_TUPOUT_ */
/* #define _PHCASEEDING_CHAIN_FORKS_ */
/* #define _PHCASEEDING_TIMER_OUT_ */

#include "ALICEKF.h"
#include "PHTrackSeeding.h"  // for PHTrackSeeding

#include <tpc/TpcGlobalPositionWrapper.h>

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
#include <memory>
#include <set>
#include <string>         // for string
#include <unordered_map>  // for map
#include <unordered_set>
#include <utility>  // for pair
#include <vector>   // for vector

#include <TFile.h>
#include <TNtuple.h>

class ActsGeometry;
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
  using pointKey = std::pair<point, TrkrDefs::cluskey>;                 // phi and z in the key
  using coordKey = std::pair<std::array<float, 2>, TrkrDefs::cluskey>;  // just use phi and Z, no longer needs the layer

  using keyList = std::vector<TrkrDefs::cluskey>;
  using keyLists = std::vector<keyList>;
  using keyListPerLayer = std::array<keyList, _NLAYERS_TPC>;
  using keySet = std::set<TrkrDefs::cluskey>;

  using keyLink = std::pair<TrkrDefs::cluskey, TrkrDefs::cluskey>;
  using keyLinks = std::vector<keyLink>;
  using keyLinkPerLayer = std::array<std::vector<keyLink>, _NLAYERS_TPC>;

  using PositionMap = std::unordered_map<TrkrDefs::cluskey, Acts::Vector3>;

  std::array<float, 55> dZ_per_layer{};
  std::array<float, 55> dphi_per_layer{};

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
      float neighbor_z_width = .01,
      float maxSinPhi = 0.999
      /* float cosTheta_limit = -0.8 */
  );

  ~PHCASeeding() override {}
  void SetSplitSeeds(bool opt = true) { _split_seeds = opt; }
  void SetLayerRange(unsigned int layer_low, unsigned int layer_up)
  {
    _start_layer = layer_low;
    _end_layer = layer_up;
  }
  void SetSearchWindow(float z_width, float phi_width)
  {
    _neighbor_z_width = z_width;
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
  {
    _min_clusters_per_seed = minClus;
    _max_clusters_per_seed = maxClus;
    if (_min_clusters_per_seed < 3)
    {
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
  void reject_zsize1_clusters(bool mode){_reject_zsize1 = mode;}
  void setNeonFraction(double frac) { Ne_frac = frac; };
  void setArgonFraction(double frac) { Ar_frac = frac; };
  void setCF4Fraction(double frac) { CF4_frac = frac; };
  void setNitrogenFraction(double frac) { N2_frac = frac; };
  void setIsobutaneFraction(double frac) { isobutane_frac = frac; };

 protected:
  int Setup(PHCompositeNode* topNode) override;
  int Process(PHCompositeNode* topNode) override;
  int InitializeGeometry(PHCompositeNode* topNode);
  /* int FindSeedsLayerSkip(double cosTheta_limit); */
  int End() override;

 private:
  bool _save_clus_proc = false;
  TFile* _f_clustering_process = nullptr;
  int _tupout_count = -1;
  int _n_tupchains = -1;
  // Save the steps of the clustering
  TNtuple* _tupclus_all = nullptr;          // all clusters
  TNtuple* _tupclus_links = nullptr;        //  clusters which are linked
  TNtuple* _tupclus_bilinks = nullptr;      // bi-linked clusters
  TNtuple* _tupclus_seeds = nullptr;        // seed (outermost three bi-linked chain
  TNtuple* _tupclus_grown_seeds = nullptr;  // seeds saved out
                                            //
  TNtuple* _tup_chainbody = nullptr;
  TNtuple* _tup_chainfork = nullptr;

  TNtuple* _tupwin_link = nullptr;       // cluster pairs with dphi and dz (very wide windows) x0,x1,y0,y1,z0,z1
  TNtuple* _tupwin_cos_angle = nullptr;  // directed cosine (all cluster triples which are found) x0,x1,x2
  TNtuple* _tupwin_seed23 = nullptr;     // from third to fourth link -- can only use passing seeds
  TNtuple* _tupwin_seedL1 = nullptr;     // from third to fourth link -- can only use passing seeds

  TNtuple* _search_windows = nullptr;  // This is really just a lazy way to store what search paramaters where used.
                                       // It would be equally valid in a TMap or map<string,float>
  // functions used to fill tuples -- only defined if _PHCASEEDING_CLUSTERLOG_TUPOUT_ is defined in preprocessor
  void write_tuples();
  void fill_tuple(TNtuple*, float, TrkrDefs::cluskey, const Acts::Vector3&) const;
  void fill_tuple_with_seed(TNtuple*, const keyList&, const PositionMap&) const;
  void process_tupout_count();
  void FillTupWinLink(bgi::rtree<pointKey, bgi::quadratic<16>>&, const coordKey&, const PositionMap&) const;
  void FillTupWinCosAngle(const TrkrDefs::cluskey, const TrkrDefs::cluskey, const TrkrDefs::cluskey, const PositionMap&, double cos_angle, bool isneg) const;
  void FillTupWinGrowSeed(const keyList& seed, const keyLink& link, const PositionMap& globalPositions) const;
  void fill_split_chains(const keyList& chain, const keyList& keylinks, const PositionMap& globalPositions, int& nchains) const;
  /* void fill_tuple_with_seed(TN */

  // have a comparator to search vector of sorted bilinks for links with starting
  // a key of a given value
  class CompKeyToBilink
  {
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
  PHCASeeding::keyLists FollowBiLinks(const keyLinks& trackSeedPairs, const keyLinkPerLayer& bilinks, const PositionMap& globalPositions) const;
  std::vector<coordKey> FillTree(bgi::rtree<pointKey, bgi::quadratic<16>>&, const keyList&, const PositionMap&, int layer);
  int FindSeedsWithMerger(const PositionMap&, const keyListPerLayer&);

  void QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>>& rtree, double phimin, double zmin, double phimax, double zmax, std::vector<pointKey>& returned_values) const;
  std::vector<TrackSeed_v2> RemoveBadClusters(const std::vector<keyList>& seeds, const PositionMap& globalPositions) const;
  double getMengerCurvature(TrkrDefs::cluskey a, TrkrDefs::cluskey b, TrkrDefs::cluskey c, const PositionMap& globalPositions) const;

  void publishSeeds(const std::vector<TrackSeed_v2>& seeds) const;

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
  unsigned int _max_clusters_per_seed = 6;  // currently not used
  unsigned int _min_clusters_per_seed = 6;  // currently the abs. number used
  //  float _cluster_z_error;
  //  float _cluster_alice_y_error;
  float _neighbor_phi_width;
  float _neighbor_z_width;
  float _clusadd_delta_dzdr_window = 0.5;
  float _clusadd_delta_dphidr2_window = 0.005;
  float _max_sin_phi;
  /* float _cosTheta_limit; */
  double _rz_outlier_threshold = 0.1;
  double _xy_outlier_threshold = 0.1;
  double _fieldDir = -1;
  bool _use_const_field = false;
  bool _split_seeds = true;
  bool _reject_zsize1 = false;
  float _const_field = 1.4;
  bool _use_fixed_clus_err = false;
  bool _pp_mode = false;
  std::array<double, 3> _fixed_clus_err = {.1, .1, .1};

  std::string m_magField;

  /// acts geometry
  ActsGeometry* m_tGeometry{nullptr};

  /// global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  std::unique_ptr<ALICEKF> fitter;

  std::unique_ptr<PHTimer> t_seed;
  std::unique_ptr<PHTimer> t_fill;
  std::unique_ptr<PHTimer> t_makebilinks;
  std::unique_ptr<PHTimer> t_makeseeds;
  /* std::array<bgi::rtree<pointKey, bgi::quadratic<16>>, _NLAYERS_TPC> _rtrees; */
  std::array<bgi::rtree<pointKey, bgi::quadratic<16>>, 3> _rtrees;  // need three layers at a time

  double Ne_frac = 0.00;
  double Ar_frac = 0.75;
  double CF4_frac = 0.20;
  double N2_frac = 0.00;
  double isobutane_frac = 0.05;
};

#endif
