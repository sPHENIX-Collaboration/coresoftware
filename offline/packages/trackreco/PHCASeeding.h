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

  // structure to hold the geometry of each cluster
  // R is calculated on the fly, as needed, and not stored here
  // phi is expensive, and used multiple times per cluster
  struct keyPoint {
    // data: the key and location geometry; R is used later (often only once) and therefore
    //       calculated on the fly as needed
    TrkrDefs::cluskey key;
    double x, y, z, phi;

    keyPoint(const TrkrDefs::cluskey _key, const Acts::Vector3& pos) : 
      key { _key }, x {pos.x()}, y{pos.y()}, z{pos.z()},
      phi{atan2(y,x)} {}; // phi range: -pi, pi
  };


  using keyPoints  = std::vector<keyPoint>;

  using keyPtr = keyPoint*;
  using keyPtrList  = std::vector<keyPtr>;
  using keyPtrLists = std::vector<std::vector<keyPtr>>;
  using keyPtrArr   = std::array<keyPtrList, _NLAYERS_TPC>;

  using keyPtrSet = std::set<keyPtr>;

  using linkIter     = std::pair<keyPtr, keyPtr>;
  using linkIterList = std::vector<linkIter>;
  using linkIterArr  = std::array<linkIterList, _NLAYERS_TPC>;

  // objects to use with boost rtree
  // the rtree pairs with pointers to the keyPoints
  using point = bg::model::point<float, 2, bg::cs::cartesian>;
  using box   = bg::model::box<point>;
  using rtreePair     = std::pair<point, keyPtr>; // phi and z in the key
  using boost_rtree   = bgi::rtree<rtreePair, bgi::quadratic<16>>;
  using rtreePairList = std::vector<rtreePair>;


  // struct to efficiency keep tracks of growing seed parameters
  // this avoid recalculating cluster values as adding clusters
  // grows outwards, always using the last three clusters
  struct Checker_dphidz
  {
   Checker_dphidz(
      const float& _delta_dzdr_window,
      const float& _delta_dphidr2_window,
      keyPtrList& seed_triplet);

    const float& delta_dzdr_window;
    const float& delta_dphidr2_window;


    std::array<float,4> z       {};
    std::array<float,4> phi     {};
    std::array<float,4> R       {};
    std::array<float,4> dR      {};
    std::array<float,4> dZdR    {};
    std::array<float,4> dphidR2 {};

    unsigned int index{0};
    int i1{1}, i2{2}, i3{3};

    inline float calc_dzdr()     { return dZdR[i3]-dZdR[i2]; };
    inline float calc_d2phidr2() { return dphidR2[i3]-2*dphidR2[i2]+dphidR2[i1]; };
    
    bool check_cluster(const keyPtr); // see if a new cluster passes dzdr and d2phidr2 
    void add_cluster(const keyPtr = nullptr); // add a cluster to end of the chain
                                              // if no entry, keep last cluster checked
    void add_clusters(const keyPtrList&);  // add the average cluster value 

    private:
    void update(const keyPtr);
    void update(const keyPtrList&);
  };

  using keyList = std::vector<TrkrDefs::cluskey>;
  using keyLists = std::vector<keyList>; 

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
  void SetMultClustersPerLayer(bool opt = true) { _doubles_in_seed = opt; }
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
  TNtuple* _tupclus_pub_seeds = nullptr;    // published seeds
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
  void FillTupWinLink(boost_rtree&, const keyPtr) const;
  void FillTupWinCosAngle(const keyPtr, const keyPtr, const keyPtr, double, bool isneg) const;
  void write_tuples();
  void fill_tuple(TNtuple*, float, keyPtr) const;
  void fill_tuple_with_seed(TNtuple*, const keyPtrList&) const;
  void process_tupout_count();
  void FillTupWinGrowSeed(const keyPtrList& seed, const keyPtr& link) const;
  void fill_split_chains(const keyPtrList& chain, const keyPtrList& keylinks, int& nchains) const;
  /* void fill_tuple_with_seed(TN */

  TpcDistortionCorrection m_distortionCorrection;

  keyPoints _cluster_pts{};
  /// get global position for a given cluster
  /**
   * uses ActsTransformation to convert cluster local position into global coordinates
   * incorporates TPC distortion correction, if present
   */

  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey, TrkrCluster*) const;

  // main segments of seeder code
  keyPtrArr GetKeyPoints(); // new
  std::pair<linkIterList, linkIterArr> CreateBilinks(keyPtrArr&);
  keyPtrLists MakeSeedTripletHeads(const linkIterList&, const linkIterArr&) const;
  void GrowSeeds(keyPtrLists&, const linkIterArr&);
  std::vector<TrackSeed_v2> RemoveBadClusters(const keyPtrLists& seeds) const;
  keyPtrList FillTree(boost_rtree&, const keyPtrList&);
  void PublishSeeds(std::vector<TrackSeed_v2>&);

  // sub functions
  void QueryTree(const boost_rtree& rtree, const keyPtr& it, const double& delta_phi, const double& delta_z, rtreePairList& returned_values) const;
  double getMengerCurvature(keyPtr a, keyPtr b, keyPtr c) const;


  // internal data
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
  bool _doubles_in_seed = false;
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
  std::unique_ptr<PHTimer> t_process;  // time steps over the major subroutines
  std::array<boost_rtree, 3> _rtrees;  // need three layers at a time

  double Ne_frac = 0.00;
  double Ar_frac = 0.75;
  double CF4_frac = 0.20;
  double N2_frac = 0.00;
  double isobutane_frac = 0.05;
};

#endif
