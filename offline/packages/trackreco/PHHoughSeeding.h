/*!
 *  \file PHHoughSeeding.h
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

#ifndef TRACKRECO_PHHOUGHSEEDING_H
#define TRACKRECO_PHHOUGHSEEDING_H

#include "PHTrackSeeding.h"

// Helix Hough includes
#ifndef __CINT__
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/VertexFinder.h>
#include <HelixHough/sPHENIXSeedFinder.h>
#endif

// PHGenFit
#include <phgenfit/Fitter.h>
//#include <phgenfit/Measurement.h>

// trackbase_historic includes
#include <trackbase_historic/SvtxTrackState.h>

#include <g4bbc/BbcVertexMap.h>

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHTimeServer.h>
#include <phool/PHTimer.h>

// standard includes
#include <float.h>
#include <list>
#include <map>
#include <memory>
#include <vector>

// forward declarations
class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class SvtxVertex;
class PHG4CellContainer;
class PHG4CylinderGeomContainer;

class PHG4HitContainer;

class TNtuple;
class TFile;

namespace PHGenFit
{
class Fitter;
class Track;
class Measurement;
} /* namespace PHGenFit */

namespace genfit
{
class GFRaveVertexFactory;
} /* namespace genfit */

/// \class PHHoughSeeding
///
/// \brief A fun4all implementation of Alan's Hough Transform
///
/// This module run Alan's Hough Transform and quick vertex finding
/// on the SvxClusterList of the event. It will produce both
/// SvxTrack and SvxSegments for the time being.
///
class PHHoughSeeding : public PHTrackSeeding
{
 public:
  PHHoughSeeding(
      const std::string &name = "PHHoughSeeding",
      unsigned int nlayers_maps = 3,
      unsigned int nlayers_intt = 8,
      unsigned int nlayers_tpc = 60,
      unsigned int seeding_nlayer = 7,
      unsigned int min_seeding_nlayer = 4);

  virtual ~PHHoughSeeding()
  {
  }

 protected:
  int Setup(PHCompositeNode *topNode);

  int Process();

  int End();

 private:
  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);

 public:
  void set_mag_field(float magField)
  {
    _magField = magField;
  }
  float get_mag_field() const
  {
    return _magField;
  }

  /// set the tracking pt-dependent chi2 cut for fast fit, cut = min(par0 + par1 / pt, max)
  void set_chi2_cut_fast(double cut_par0, double cut_par1 = 0.0,
                         double cut_max = FLT_MAX)
  {
    _chi2_cut_fast_par0 = cut_par0;
    _chi2_cut_fast_par1 = cut_par1;
    _chi2_cut_fast_max = cut_max;
  }

  /// set the tracking chi2 cut for full fit
  void set_chi2_cut_full(double chi2_cut)
  {
    _chi2_cut_full = chi2_cut;
  }
  /// get the tracking chi2 cut for full fit
  double get_chi2_cut_full() const
  {
    return _chi2_cut_full;
  }

  /// set early combination-land chi2 cut(?)
  void set_ca_chi2_cut(double chi2_cut)
  {
    _ca_chi2_cut = chi2_cut;
  }
  /// get early combination-land chi2 cut(?)
  double get_ca_chi2_cut() const
  {
    return _ca_chi2_cut;
  }

  /// set early curvature cut between hits, lower values are more open
  void set_cos_angle_cut(double cos_angle_cut)
  {
    _cos_angle_cut = cos_angle_cut;
  }
  /// get early curvature cut between hits, lower values are more open
  double get_cos_angle_cut() const
  {
    return _cos_angle_cut;
  }

  void set_min_pT(float pt)
  {
    _min_pt = pt;
  }

  /// set the z0 search window
  void set_z0_range(float min_z0, float max_z0)
  {
    _min_z0 = min_z0;
    _max_z0 = max_z0;
  }

  /// set the r search window
  void set_r_max(float max_r)
  {
    _max_r = max_r;
  }

  /// radiation length per layer, sequential layer indexes required here
  void set_material(int layer, float value);

  /// set internal ghost rejection
  void setRejectGhosts(bool rg)
  {
    _reject_ghosts = rg;
  }
  /// set for internal hit rejection
  void setRemoveHits(bool rh)
  {
    _remove_hits = rh;
  }

  /// adjusts the rate of zooming
  void setBinScale(float scale)
  {
    _bin_scale = scale;
  }
  /// adjusts the rate of zooming
  void setZBinScale(float scale)
  {
    _z_bin_scale = scale;
  }

  /// turn on DCA limitation
  void setCutOnDCA(bool cod)
  {
    _cut_on_dca = cod;
  }
  /// sets an upper limit on X-Y DCA
  void setDCACut(float dcut)
  {
    _dcaxy_cut = dcut;
  }
  /// sets an upper limit on Z DCA
  void setDCAZCut(float dzcut)
  {
    _dcaz_cut = dzcut;
  }

  /// adjust the fit pt by a recalibration factor (constant B versus real mag
  /// field)
  void setPtRescaleFactor(float pt_rescale)
  {
    _pt_rescale = pt_rescale;
  }

  /// adjust the relative voting error scale w.r.t. the cluster size
  void setVoteErrorScale(unsigned int layer, float scale)
  {
    if (scale > 0.0)
    {
      _vote_error_scale.at(layer) = scale;
    }
    else
    {
      std::cout << "PHHoughSeeding::setVoteErrorScale : scale must be "
                   "greater than zero ... doing nothing"
                << std::endl;
    }
  }
  /// adjust the relative fit error scale w.r.t. the cluster size
  void setFitErrorScale(unsigned int layer, float scale)
  {
    if (scale > 0.0)
    {
      _fit_error_scale.at(layer) = scale;
    }
    else
    {
      std::cout << "PHHoughSeeding::setFitErrorScale : scale must be "
                   "greater than zero ... doing nothing"
                << std::endl;
    }
  }

  //---deprecated---------------------------------------------------------------

  /// set option to produce initial vertex for further tracking
  void set_use_vertex(bool b)
  {
  }

  void setInitialResMultiplier(int beta)
  {
  }
  void setFullResMultiplier(int lambda)
  {
  }

  /// set the minimum pT to try to find during initial vertex finding tracking
  void set_min_pT_init(float PT)
  {
  }

  /// limit the maximum error reported by cluster (in number of cell units)
  void setMaxClusterError(float max_cluster_error)
  {
  }

  /// use the cell size as cluster size instead of value stored on cluster
  void setUseCellSize(bool use_cell_size)
  {
  }

  /// set the tracking chi2 for initial vertex finding
  void set_chi2_cut_init(double chi2_cut)
  {
  }

  const std::vector<int> &get_seeding_layer() const
  {
    return _seeding_layer;
  }

  //	void set_seeding_layer(const std::vector<int>& seedingLayer) {
  //		_seeding_layer = seedingLayer;
  //	}

  void set_seeding_layer(const int *seedingLayer, const int n)
  {
    _seeding_layer.clear();
    _seeding_layer.assign(seedingLayer, seedingLayer + n);
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

  void set_vertex_error(const float a)
  {
    _vertex_error.clear();
    _vertex_error.assign(3, a);
  }

  unsigned int get_min_nlayers_seeding() const
  {
    return _min_nlayers_seeding;
  }

  void set_min_nlayers_seeding(unsigned int minNlayersSeeding)
  {
    _min_nlayers_seeding = minNlayersSeeding;
    _min_combo_hits = minNlayersSeeding;
  }

#ifndef __CINT__

 private:
  //--------------
  // InitRun Calls
  //--------------

  /// create new node output pointers
  //int CreateNodes(PHCompositeNode *topNode);

  /// scan tracker geometry objects
  int InitializeGeometry(PHCompositeNode *topNode);

  /// track propagation
  int InitializePHGenFit(PHCompositeNode *topNode);

  /// code to setup seed tracking objects
  int setup_seed_tracker_objects();

  /// code to setup initial vertexing tracker
  int setup_initial_tracker_object();

  /// code to setup full tracking object
  int setup_tracker_object();

  //--------------------
  // Process Event Calls
  //--------------------

  /// fetch node pointers
  //int GetNodes(PHCompositeNode *topNode);

  /// code to translate into the HelixHough universe
  int translate_input();

  /// code to combine seed tracking vertex with BBCZ if available
  int fast_composite_seed();

  /// code to seed vertex from bbc
  int fast_vertex_from_bbc();

  /// code to seed vertex from initial tracking using a broad search window
  int fast_vertex_guessing();

  /// code to produce an initial track vertex from the seed
  int initial_vertex_finding();

  /// SvtxVtxMap[0] -> _vertex and _vertex_error
  int vertexing();

  /// code to perform the final tracking and vertexing
  int full_track_seeding();

  /// code to translate back to the SVTX universe
  int export_output();

  //!
  int CleanupSeedsByHitPattern();

  //!
  int CleanupTracksByHitPattern();

  //!
  int CleanupSeeds();

  //!
  void print_timers();

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
                                                   Eigen::Matrix<float, 5, 5> const &input,
                                                   Eigen::Matrix<float, 6, 6> &output);

  /// translate the clusters, tracks, and vertex from one origin to another
  void shift_coordinate_system(double dx, double dy, double dz);

  int _event;
  PHTimer *_t_seeding;
  PHTimer *_t_seed_init1;
  PHTimer *_t_seed_init2;
  PHTimer *_t_seed_init3;
  PHTimer *_t_seeds_cleanup;
  PHTimer *_t_translate_to_PHGenFitTrack;
  PHTimer *_t_translate1;
  PHTimer *_t_translate2;
  PHTimer *_t_translate3;
  PHTimer *_t_kalman_pat_rec;
  PHTimer *_t_search_clusters;
  PHTimer *_t_search_clusters_encoding;
  PHTimer *_t_search_clusters_map_iter;
  PHTimer *_t_track_propagation;
  PHTimer *_t_full_fitting;
  PHTimer *_t_output_io;

  std::vector<int> _seeding_layer;  //layer numbers that are used for seeding

  unsigned int _nlayers_seeding;        ///< number of detector layers
  unsigned int _min_nlayers_seeding;    ///< minimum number of layers to make a track
  std::vector<float> _radii;            ///< radial distance of each layer (cm)
  std::vector<float> _material;         ///< material at each layer in rad. lengths
  std::map<int, float> _user_material;  ///< material in user ladder indexes

  float _magField;  ///< in Tesla

  bool _reject_ghosts;
  bool _remove_hits;

  float _min_pt;
  float _min_z0;
  float _max_z0;
  float _max_r;

  bool _cut_on_dca;
  float _dcaxy_cut;
  float _dcaz_cut;

  double _chi2_cut_fast_par0;  ///< fit quality chisq/dof for fast fit track fitting
  double _chi2_cut_fast_par1;  ///< fit quality chisq/dof for fast fit track fitting
  double _chi2_cut_fast_max;   ///< fit quality chisq/dof for fast fit track fitting
  double _chi2_cut_full;       ///< fit quality chisq/dof for kalman track fitting
  double _ca_chi2_cut;         ///< initial combination cut?
  double _cos_angle_cut;       ///< curvature restriction on cluster combos

  float _bin_scale;
  float _z_bin_scale;

  unsigned int _min_combo_hits;  ///< minimum hits to enter combination gun
  unsigned int _max_combo_hits;  ///< maximum hits to enter combination gun

  float _pt_rescale;
  std::vector<float> _fit_error_scale;
  std::vector<float> _vote_error_scale;

  /// recorded layer indexes to internal sequential indexes
  std::map<int, unsigned int> _layer_ilayer_map;

  // object storage
  std::vector<SimpleHit3D> _clusters;  ///< working array of clusters

  std::vector<SimpleTrack3D> _tracks;                      ///< working array of tracks
  std::vector<double> _track_errors;                       ///< working array of track chisq
  std::vector<Eigen::Matrix<float, 5, 5> > _track_covars;  ///< working array of track covariances

  std::vector<SimpleTrack3D> _all_tracks;                      ///< working array of tracks
  std::vector<double> _all_track_errors;                       ///< working array of track chisq
  std::vector<Eigen::Matrix<float, 5, 5> > _all_track_covars;  ///< working array of track covariances

  std::vector<float> _vertex;        ///< working array for collision vertex
  std::vector<float> _vertex_error;  ///< sqrt(cov)

  // track finding routines
  sPHENIXSeedFinder *_tracker;            ///< finds full tracks
  sPHENIXSeedFinder *_tracker_vertex;     ///< finds a subset of tracks for initial vertex-finding
  sPHENIXSeedFinder *_tracker_etap_seed;  ///< finds a subset of tracks for the vertex guess
  sPHENIXSeedFinder *_tracker_etam_seed;  ///< finds a subset of tracks for the vertex guess
  VertexFinder _vertexFinder;             ///< vertex finding object

  // node pointers
  BbcVertexMap *_bbc_vertexes;
  //SvtxTrackMap* _track_map;
  //SvtxVertexMap* _vertex_map;

  //nodes to get norm vector
  //SvtxHitMap *_svtxhitsmap;

  int *_hit_used_map;
  int _hit_used_map_size;

  PHG4CylinderGeomContainer *_geom_container_intt;
  PHG4CylinderGeomContainer *_geom_container_maps;

  bool _analyzing_mode;
  TFile *_analyzing_file;
  TNtuple *_analyzing_ntuple;

  //! Cleanup Seeds
  float _max_merging_dphi;
  float _max_merging_deta;
  float _max_merging_dr;
  float _max_merging_dz;

  //! if two seeds have common hits more than this number, merge them
  unsigned int _max_share_hits;

  double _cut_min_pT;

  unsigned int _nlayers_maps;
  unsigned int _nlayers_intt;
  unsigned int _nlayers_tpc;

  int _nlayers_all;

  std::map<int, unsigned int> _layer_ilayer_map_all;
  std::vector<float> _radii_all;

#endif  // __CINT__
};

#endif
