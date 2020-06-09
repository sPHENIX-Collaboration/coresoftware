/*!
 *  \file PHGenFitTrkProp.h
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

#ifndef TRACKRECO_PHGENFITTRKPROP_H
#define TRACKRECO_PHGENFITTRKPROP_H

#include "PHTrackPropagating.h"

#include <trackbase_historic/SvtxTrackMap.h>

#include <trackbase/TrkrDefs.h>               // for cluskey

#include <Eigen/Core>                         // for Matrix
// needed, it crashes on Ubuntu using singularity with local cvmfs install
// shared pointer later on uses this, forward declaration does not cut it
#include <phgenfit/Track.h> 
#include <gsl/gsl_rng.h>

// standard includes
#include <array>
#include <list>
#include <map>
#include <memory>
#include <ostream>                            // for basic_ostream::operator<<
#include <string>                             // for string
#include <utility>                            // for pair
#include <vector>

// forward declarations
class BbcVertexMap;
class PHCompositeNode;
class PHG4CylinderGeomContainer;
class PHTimer;
class SvtxTrack;
class TrkrCluster;

class TFile;
class TNtuple;

namespace PHGenFit
{
class Fitter;
class Measurement;
} /* namespace PHGenFit */

///
/// \class PHGenFitTrkProp
/// \brief Propagate tracklet to full track using GenFit
///
class PHGenFitTrkProp : public PHTrackPropagating
{
 public:
  PHGenFitTrkProp(
      const std::string& name = "PHGenFitTrkProp",
      unsigned int nlayers_maps = 3,
      unsigned int nlayers_intt = 8,
      unsigned int nlayers_tpc = 60,
      unsigned int nlayers_micromegas = 0);

 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process() override;

  int End() override;

 private:
  /// fetch node pointers
  int GetNodes(PHCompositeNode* topNode);

 public:
  //	int Init(PHCompositeNode *topNode);
  //	int InitRun(PHCompositeNode *topNode);
  //	int process_event(PHCompositeNode *topNode);
  //	int End(PHCompositeNode *topNode);

  struct TrackQuality
  {
    int nhits = 0;
    float chi2 = 0;
    int ndf = 0;
    int nmicromegas = 0;
    int ntpc = 0;
    int nintt = 0;
    int nmaps = 0;

    TrackQuality(int nhits_, float chi2_, int ndf_)
      : nhits(nhits_)
      , chi2(chi2_)
      , ndf(ndf_)
    {
    }

    TrackQuality(int nhits_, float chi2_, int ndf_, int nmicromegas_, int ntpc_, int nintt_, int nmaps_)
      : nhits(nhits_)
      , chi2(chi2_)
      , ndf(ndf_)
      , nmicromegas(nmicromegas_)
      , ntpc(ntpc_)
      , nintt(nintt_)
      , nmaps(nmaps_)
    {
    }

    bool operator<(const TrackQuality& b) const
    {
      if (nhits != b.nhits)
        return nhits > b.nhits;
      else
        return chi2 / ndf < b.chi2 / b.ndf;
    }

    friend std::ostream& operator<<(std::ostream& os, const PHGenFitTrkProp::TrackQuality& tq)
    {
      os
          << tq.nhits << ", "
          << tq.chi2 << ", " << tq.ndf << ", "
          << tq.nmicromegas << ", " << tq.ntpc << ", " << tq.nintt << ", " << tq.nmaps
          << std::endl;

      return os;
    }
  };

  typedef std::list<std::pair<TrackQuality, std::shared_ptr<PHGenFit::Track> > > MapPHGenFitTrack;

  float get_search_win_phi() const
  {
    return _search_win_phi;
  }

  void set_search_win_phi(float searchWinMultiplier)
  {
    _search_win_phi = searchWinMultiplier;
  }

  bool is_do_evt_display() const
  {
    return _do_evt_display;
  }

  void set_do_evt_display(bool doEvtDisplay)
  {
    _do_evt_display = doEvtDisplay;
  }

  const std::string& get_track_fitting_alg_name() const
  {
    return _track_fitting_alg_name;
  }

  void set_track_fitting_alg_name(const std::string& trackFittingAlgName)
  {
    _track_fitting_alg_name = trackFittingAlgName;
  }

  void set_analyzing_mode(bool analyzingMode)
  {
    _analyzing_mode = analyzingMode;
  }

  float get_max_merging_deta() const
  {
    return _max_merging_deta;
  }

  void set_max_merging_deta(float maxMergingDeta)
  {
    _max_merging_deta = maxMergingDeta;
  }

  float get_max_merging_dphi() const
  {
    return _max_merging_dphi;
  }

  void set_max_merging_dphi(float maxMergingDphi)
  {
    _max_merging_dphi = maxMergingDphi;
  }

  float get_max_merging_dr() const
  {
    return _max_merging_dr;
  }

  void set_max_merging_dr(float maxMergingDr)
  {
    _max_merging_dr = maxMergingDr;
  }

  float get_max_merging_dz() const
  {
    return _max_merging_dz;
  }

  void set_max_merging_dz(float maxMergingDz)
  {
    _max_merging_dz = maxMergingDz;
  }

  float get_search_win_theta() const
  {
    return _search_win_theta;
  }

  void set_search_win_theta(float searchWinZ)
  {
    _search_win_theta = searchWinZ;
  }

  float get_max_incr_chi2() const
  {
    return _max_incr_chi2;
  }

  void set_max_incr_chi2(float maxIncrChi2)
  {
    _max_incr_chi2 = maxIncrChi2;
  }

  unsigned int get_max_consecutive_missing_layer() const
  {
    return _max_consecutive_missing_layer;
  }

  void set_max_consecutive_missing_layer(
      unsigned int maxConsecutiveMissingLayer)
  {
    _max_consecutive_missing_layer = maxConsecutiveMissingLayer;
  }

  unsigned int get_min_good_track_hits() const
  {
    return _min_good_track_hits;
  }

  void set_min_good_track_hits(unsigned int minGoodTrackHits)
  {
    _min_good_track_hits = minGoodTrackHits;
  }

  unsigned int get_max_share_hits() const
  {
    return _max_share_hits;
  }

  void set_max_share_hits(unsigned int maxShareHits)
  {
    _max_share_hits = maxShareHits;
  }

  float get_max_splitting_chi2() const
  {
    return _max_splitting_chi2;
  }

  void set_max_splitting_chi2(float maxSplittingChi2)
  {
    _max_splitting_chi2 = maxSplittingChi2;
  }

  int get_init_direction() const
  {
    return _init_direction;
  }

  void set_init_direction(int initDirection)
  {
    _init_direction = initDirection;
  }

  float get_max_search_win_phi_tpc() const
  {
    return _max_search_win_phi_tpc;
  }

  void set_max_search_win_phi_tpc(float maxSearchWinPhi)
  {
    _max_search_win_phi_tpc = maxSearchWinPhi;
  }

  float get_max_search_win_theta_tpc() const
  {
    return _max_search_win_theta_tpc;
  }

  void set_max_search_win_theta_tpc(float maxSearchWinZ)
  {
    _max_search_win_theta_tpc = maxSearchWinZ;
  }
  float get_max_search_win_phi_micromegas(int layer) const
  {
    return _max_search_win_phi_micromegas[layer];
  }

  void set_max_search_win_phi_micromegas(int layer, float maxSearchWinPhi)
  {
    _max_search_win_phi_micromegas[layer] = maxSearchWinPhi;
  }

  float get_max_search_win_theta_micromegas(int layer) const
  {
    return _max_search_win_theta_micromegas[layer];
  }

  void set_max_search_win_theta_micromegas(int layer, float maxSearchWinZ)
  {
    _max_search_win_theta_micromegas[layer] = maxSearchWinZ;
  }

  float get_blowup_factor() const
  {
    return _blowup_factor;
  }

  void set_blowup_factor(float blowupFactor)
  {
    _blowup_factor = blowupFactor;
  }

  float get_max_search_win_phi_intt(int inttlayer) const
  {
    return _max_search_win_phi_intt[inttlayer];
  }

  void set_max_search_win_phi_intt(int inttlayer, float maxSearchWinPhiIntt)
  {
    _max_search_win_phi_intt[inttlayer] = maxSearchWinPhiIntt;
  }

  float get_max_search_win_phi_maps() const
  {
    return _max_search_win_phi_maps;
  }

  void set_max_search_win_phi_maps(float maxSearchWinPhiMaps)
  {
    _max_search_win_phi_maps = maxSearchWinPhiMaps;
  }

  float get_max_search_win_theta_intt(int inttlayer) const
  {
    return _max_search_win_theta_intt[inttlayer];
  }

  void set_max_search_win_theta_intt(int inttlayer, float maxSearchWinThetaIntt)
  {
    _max_search_win_theta_intt[inttlayer] = maxSearchWinThetaIntt;
  }

  float get_max_search_win_theta_maps() const
  {
    return _max_search_win_theta_maps;
  }

  void set_max_search_win_theta_maps(float maxSearchWinThetaMaps)
  {
    _max_search_win_theta_maps = maxSearchWinThetaMaps;
  }

  float get_min_search_win_phi_intt(int inttlayer) const
  {
    return _min_search_win_phi_intt[inttlayer];
  }

  void set_min_search_win_phi_intt(int inttlayer, float minSearchWinPhiIntt)
  {
    _min_search_win_phi_intt[inttlayer] = minSearchWinPhiIntt;
  }

  float get_min_search_win_phi_maps() const
  {
    return _min_search_win_phi_maps;
  }

  void set_min_search_win_phi_maps(float minSearchWinPhiMaps)
  {
    _min_search_win_phi_maps = minSearchWinPhiMaps;
  }

  float get_min_search_win_phi_tpc() const
  {
    return _min_search_win_phi_tpc;
  }

  void set_min_search_win_phi_tpc(float minSearchWinPhiTpc)
  {
    _min_search_win_phi_tpc = minSearchWinPhiTpc;
  }

  float get_min_search_win_phi_micromegas(int layer) const
  {
    return _min_search_win_phi_micromegas[layer];
  }

  void set_min_search_win_phi_micromegas(int layer, float minSearchWinPhimicromegas)
  {
    _min_search_win_phi_micromegas[layer] = minSearchWinPhimicromegas;
  }

  float get_min_search_win_theta_intt(int inttlayer) const
  {
    return _min_search_win_theta_intt[inttlayer];
  }

  void set_min_search_win_theta_intt(int inttlayer, float minSearchWinThetaIntt)
  {
    _min_search_win_theta_intt[inttlayer] = minSearchWinThetaIntt;
  }

  float get_min_search_win_theta_maps() const
  {
    return _min_search_win_theta_maps;
  }

  void set_min_search_win_theta_maps(float minSearchWinThetaMaps)
  {
    _min_search_win_theta_maps = minSearchWinThetaMaps;
  }

  float get_min_search_win_theta_tpc() const
  {
    return _min_search_win_theta_tpc;
  }

  void set_min_search_win_theta_tpc(float minSearchWinThetaTpc)
  {
    _min_search_win_theta_tpc = minSearchWinThetaTpc;
  }

  float get_min_search_win_theta_micromegas(int layer) const
  {
    return _min_search_win_theta_micromegas[layer];
  }

  void set_min_search_win_theta_micromegas(int layer, float minSearchWinThetamicromegas)
  {
    _min_search_win_theta_micromegas[layer] = minSearchWinThetamicromegas;
  }

  int get_primary_pid_guess() const
  {
    return _primary_pid_guess;
  }

  void set_primary_pid_guess(int primaryPidGuess)
  {
    _primary_pid_guess = primaryPidGuess;
  }

 private:

  //*@name utility methods
  //@{
  inline bool is_maps_layer( unsigned int layer ) const
  { return layer >= _firstlayer_maps && layer < _firstlayer_maps + _nlayers_maps; }

  inline bool is_intt_layer( unsigned int layer ) const
  { return layer >= _firstlayer_intt && layer < _firstlayer_intt + _nlayers_intt; }

  inline bool is_tpc_layer( unsigned int layer ) const
  { return layer >= _firstlayer_tpc && layer < _firstlayer_tpc + _nlayers_tpc; }

  inline bool is_micromegas_layer( unsigned int layer ) const
  { return layer >= _firstlayer_micromegas && layer < _firstlayer_micromegas + _nlayers_micromegas; }
  //@}

  //--------------
  // InitRun Calls
  //--------------

  /// Init projection r
  int InitializeGeometry(PHCompositeNode* topNode);

  /// track propagation
  int InitializePHGenFit(PHCompositeNode* topNode);

  //--------------------
  // Process Event Calls
  //--------------------

  //!
  int check_track_exists(MapPHGenFitTrack::iterator, SvtxTrackMap::Iter);

  //! Main function
  int KalmanTrkProp();

  //!
  int ExportOutput();

  //!
  void print_timers();

  //--------------------
  //
  //--------------------

  //! KalmanTrkProp Call.
  int BuildLayerZPhiHitMap(const unsigned int ivert);

  //! layer: 7 bits, z: 11 bits, phi: 14 bits
  unsigned int encode_cluster_index(const unsigned int layer, const float z, const float rphi);

  unsigned int encode_cluster_index(const unsigned int layer, const unsigned int iz, const unsigned int irphi);

  //! KalmanTrkProp Call.
  int SvtxTrackToPHGenFitTracks(const SvtxTrack* svtxtrack);

  //	int TrackPropPatRec(PHCompositeNode* topNode,
  //			//const int iPHGenFitTrack, std::shared_ptr<PHGenFit::Track> &track,
  //			MapPHGenFitTrack::iterator &track_iter,
  //			const unsigned int init_layer = 0, const unsigned int end_layer = 66,
  //			const bool use_fitted_state_once = false);
  int TrackPropPatRec(
		      const unsigned int ivert,
		      //const int iPHGenFitTrack, std::shared_ptr<PHGenFit::Track> &track,
		      MapPHGenFitTrack::iterator& track_iter,
		      const unsigned int init_layer = 0, const unsigned int end_layer = 66,
		      const bool use_fitted_state_once = false);
  
  //!
  //PHGenFit::Measurement* SvtxClusterToPHGenFitMeasurement(const SvtxCluster* cluster);
  PHGenFit::Measurement* TrkrClusterToPHGenFitMeasurement(const TrkrCluster* cluster);

  //! TrackPropPatRec Call.
  std::vector<TrkrDefs::cluskey> SearchHitsNearBy(const unsigned int ivert, const unsigned int layer, const float z_center, const float phi_center, const float z_window, const float phi_window);

  //! ExportOutput Call. Make SvtxTrack from PHGenFit::Track and set of clusters
  //std::shared_ptr<SvtxTrack> MakeSvtxTrack(const int genfit_track_ID, const SvtxVertex * vertex = NULL);
  int OutputPHGenFitTrack(MapPHGenFitTrack::iterator, SvtxTrackMap::Iter);

  //------------------
  // Subfunction Calls
  //------------------

  /// convert from inverse curvature to momentum
  float kappaToPt(float kappa);

  /// convert from momentum to inverse curvature
  float ptToKappa(float pt);

  /// convert the covariance from HelixHough coords to x,y,z,px,py,pz
  void convertHelixCovarianceToEuclideanCovariance(float B, float phi,
                                                   float d, float kappa, float z0, float dzdl,
                                                   Eigen::Matrix<float, 5, 5> const& input,
                                                   Eigen::Matrix<float, 6, 6>& output);

  /// translate the clusters, tracks, and vertex from one origin to another
  void shift_coordinate_system(double dx, double dy, double dz);

  int _event = 0;
  PHTimer* _t_seeds_cleanup = nullptr;
  PHTimer* _t_translate_to_PHGenFitTrack = nullptr;
  PHTimer* _t_translate1 = nullptr;
  PHTimer* _t_translate2 = nullptr;
  PHTimer* _t_translate3 = nullptr;
  PHTimer* _t_kalman_pat_rec = nullptr;
  PHTimer* _t_search_clusters = nullptr;
  PHTimer* _t_search_clusters_encoding = nullptr;
  PHTimer* _t_search_clusters_map_iter = nullptr;
  PHTimer* _t_track_propagation = nullptr;
  PHTimer* _t_full_fitting = nullptr;
  PHTimer* _t_output_io = nullptr;
  
  // object storage

  std::vector<std::vector<float>> _vertex;        ///< working array for collision vertex list
  std::vector<std::vector<float>> _vertex_error;  ///< sqrt(cov)
    
  // node pointers
  BbcVertexMap* _bbc_vertexes = nullptr;

  //nodes to get norm vector
  //SvtxHitMap* _svtxhitsmap;

  int* _hit_used_map = nullptr;
  int _hit_used_map_size = 0;
  std::multimap<TrkrDefs::cluskey, unsigned int> _gftrk_hitkey_map;

  PHG4CylinderGeomContainer* _geom_container_intt = nullptr;
  PHG4CylinderGeomContainer* _geom_container_maps = nullptr;

  bool _analyzing_mode = false;
  TFile* _analyzing_file = nullptr;
  TNtuple* _analyzing_ntuple = nullptr;

  //! Cleanup Seeds
  float _max_merging_dphi = 0.1;
  float _max_merging_deta = 0.1;
  float _max_merging_dr = 0.1;
  float _max_merging_dz = 0.1;

  //! if two seeds have common hits more than this number, merge them
  unsigned int _max_share_hits = 3;

  std::unique_ptr<PHGenFit::Fitter> _fitter;

  //! KalmanFitterRefTrack, KalmanFitter, DafSimple, DafRef
  std::string _track_fitting_alg_name = "KalmanFitter";

  int _primary_pid_guess = 211;
  double _cut_min_pT = 0.2;

  bool _do_evt_display = false;

  unsigned int _nlayers_maps = 3;
  unsigned int _nlayers_intt = 4;
  unsigned int _nlayers_tpc = 48;
  unsigned int _nlayers_micromegas = 0;

  int _nlayers_all = 55;

  unsigned int _firstlayer_maps = 0;
  unsigned int _firstlayer_intt = 3;
  unsigned int _firstlayer_tpc = 7;
  unsigned int _firstlayer_micromegas = 55;

  std::map<int, unsigned int> _layer_ilayer_map_all;
  std::vector<float> _radii_all;

  // TODO: might need to use layer dependent windows because micromegas are 1D measurements
  std::array<float,2> _max_search_win_phi_micromegas = {{ 0.004, 0.62}};
  std::array<float,2> _min_search_win_phi_micromegas = {{ 0, 0.31 }};
  std::array<float,2> _max_search_win_theta_micromegas = {{ 1.24, 0.004}};
  std::array<float,2> _min_search_win_theta_micromegas = {{ 0.62, 0}};

  float _max_search_win_phi_tpc = 0.004;
  float _min_search_win_phi_tpc = 0;
  float _max_search_win_theta_tpc = 0.004;
  float _min_search_win_theta_tpc = 0;

  std::array<float,4> _max_search_win_phi_intt = {{ 0.2, 0.2, 0.005, 0.005 }};
  std::array<float,4> _min_search_win_phi_intt = {{ 0.2, 0.2, 0, 0 }};
  std::array<float,4> _max_search_win_theta_intt = {{ 0.01, 0.01, 0.2, 0.2 }};
  std::array<float,4> _min_search_win_theta_intt = {{ 0, 0, 0.2, 0.2 }};

  float _max_search_win_phi_maps = 0.005;
  float _min_search_win_phi_maps = 0;
  float _max_search_win_theta_maps = 0.04;
  float _min_search_win_theta_maps = 0;

  float _search_win_phi = 20;
  float _search_win_theta = 20;
  std::map<int, float> _search_wins_phi;
  std::map<int, float> _search_wins_theta;

  std::vector<std::multimap<unsigned int, TrkrDefs::cluskey>> _layer_thetaID_phiID_clusterID;

  float _half_max_theta = M_PI/2;
  float _half_max_phi = M_PI;
  float _layer_thetaID_phiID_clusterID_phiSize = 0.12/30;
  float _layer_thetaID_phiID_clusterID_zSize = 0.17/30;

  MapPHGenFitTrack _PHGenFitTracks;
  //! +1: inside out; -1: outside in
  int _init_direction = -1;
  float _blowup_factor = 1;

  unsigned int _max_consecutive_missing_layer = 20;
  float _max_incr_chi2 = 20;
  std::map<int, float> _max_incr_chi2s;

  float _max_splitting_chi2 = 20;

  unsigned int _min_good_track_hits = 30;

};

#endif
