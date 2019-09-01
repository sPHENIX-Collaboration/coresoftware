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

#if !defined(__CINT__) || defined(__CLING__)
#include <Eigen/Core>                         // for Matrix
#endif

// standard includes
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
class Track;
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
      unsigned int nlayers_tpc = 60);

  virtual ~PHGenFitTrkProp();

 protected:
  int Setup(PHCompositeNode* topNode);

  int Process();

  int End();

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
    int nhits;
    float chi2;
    int ndf;

    int ntpc;
    int nintt;
    int nmaps;

    TrackQuality(int nhits_, float chi2_, int ndf_)
      : nhits(nhits_)
      , chi2(chi2_)
      , ndf(ndf_)
      , ntpc(0)
      , nintt(0)
      , nmaps(0)
    {
    }

    TrackQuality(int nhits_, float chi2_, int ndf_, int ntpc_, int nintt_, int nmaps_)
      : nhits(nhits_)
      , chi2(chi2_)
      , ndf(ndf_)
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
          << tq.ntpc << ", " << tq.nintt << ", " << tq.nmaps
          << std::endl;

      return os;
    }
  };

#ifndef __CINT__
  typedef std::list<std::pair<TrackQuality, std::shared_ptr<PHGenFit::Track> > > MapPHGenFitTrack;
#endif

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

  void set_vertex_error(const float a)
  {
    _vertex_error.clear();
    // _vertex_error.assign(3, a);
  }

  int get_primary_pid_guess() const
  {
    return _primary_pid_guess;
  }

  void set_primary_pid_guess(int primaryPidGuess)
  {
    _primary_pid_guess = primaryPidGuess;
  }

#ifndef __CINT__

 private:
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

  int _event;
  PHTimer* _t_seeds_cleanup;
  PHTimer* _t_translate_to_PHGenFitTrack;
  PHTimer* _t_translate1;
  PHTimer* _t_translate2;
  PHTimer* _t_translate3;
  PHTimer* _t_kalman_pat_rec;
  PHTimer* _t_search_clusters;
  PHTimer* _t_search_clusters_encoding;
  PHTimer* _t_search_clusters_map_iter;
  PHTimer* _t_track_propagation;
  PHTimer* _t_full_fitting;
  PHTimer* _t_output_io;

  // object storage

  std::vector<std::vector<float>> _vertex;        ///< working array for collision vertex list
  std::vector<std::vector<float>> _vertex_error;  ///< sqrt(cov)
  
  //  std::vector<float> _vertex;        ///< working array for collision vertex
  // std::vector<float> _vertex_error;  ///< sqrt(cov)
  
  // node pointers
  BbcVertexMap* _bbc_vertexes;

  //nodes to get norm vector
  //SvtxHitMap* _svtxhitsmap;

  int* _hit_used_map;
  int _hit_used_map_size;
  std::multimap<TrkrDefs::cluskey, unsigned int> _gftrk_hitkey_map;
  // PHG4CellContainer* _cells_svtx;
  // PHG4CellContainer* _cells_intt;
  //PHG4CellContainer* _cells_maps;

  PHG4CylinderGeomContainer* _geom_container_intt;
  PHG4CylinderGeomContainer* _geom_container_maps;

  bool _analyzing_mode;
  TFile* _analyzing_file;
  TNtuple* _analyzing_ntuple;

  //! Cleanup Seeds
  float _max_merging_dphi;
  float _max_merging_deta;
  float _max_merging_dr;
  float _max_merging_dz;

  //! if two seeds have common hits more than this number, merge them
  unsigned int _max_share_hits;

  PHGenFit::Fitter* _fitter;

  //! KalmanFitterRefTrack, KalmanFitter, DafSimple, DafRef
  //PHGenFit::Fitter::FitterType _track_fitting_alg_name;
  std::string _track_fitting_alg_name;

  int _primary_pid_guess;
  double _cut_min_pT;

  bool _do_evt_display;

  unsigned int _nlayers_maps;
  unsigned int _nlayers_intt;
  unsigned int _nlayers_tpc;

  int _nlayers_all;

  std::map<int, unsigned int> _layer_ilayer_map_all;
  std::vector<float> _radii_all;

  float _max_search_win_phi_tpc;
  float _min_search_win_phi_tpc;
  float _max_search_win_theta_tpc;
  float _min_search_win_theta_tpc;

  float _max_search_win_phi_intt[8];
  float _min_search_win_phi_intt[8];
  float _max_search_win_theta_intt[8];
  float _min_search_win_theta_intt[8];

  float _max_search_win_phi_maps;
  float _min_search_win_phi_maps;
  float _max_search_win_theta_maps;
  float _min_search_win_theta_maps;

  float _search_win_phi;
  float _search_win_theta;
  std::map<int, float> _search_wins_phi;
  std::map<int, float> _search_wins_theta;

  std::vector<std::multimap<unsigned int, TrkrDefs::cluskey>> _layer_thetaID_phiID_clusterID;

  float _half_max_theta;
  float _half_max_phi;
  float _layer_thetaID_phiID_clusterID_phiSize;
  float _layer_thetaID_phiID_clusterID_zSize;

  MapPHGenFitTrack _PHGenFitTracks;
  //! +1: inside out; -1: outside in
  int _init_direction;
  float _blowup_factor;

  unsigned int _max_consecutive_missing_layer;
  float _max_incr_chi2;
  std::map<int, float> _max_incr_chi2s;

  float _max_splitting_chi2;

  unsigned int _min_good_track_hits;

#endif  // __CINT__
};

#endif
