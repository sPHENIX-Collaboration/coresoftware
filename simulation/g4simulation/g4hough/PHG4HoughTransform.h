#ifndef __PHG4HOUGHTRANSFORM_H__
#define __PHG4HOUGHTRANSFORM_H__

//===========================================================
/// \file PHG4HoughTransform.h
/// \brief A fun4all implementation of Alan's Hough Transform
/// \author Matt Wysocki (copied from SvxHoughTransform)
/// go to
/// https://www.phenix.bnl.gov/WWW/offline/wikioffline/index.php/SvxHoughTransform
/// \edited by Theo Koblesky to conform to new changes in HelixHough (1/10/2012)

// edited week of 11/14/2012 by Alan Dion to use new HelixHough
// classes

//===========================================================

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHTimeServer.h>
#include <g4bbc/BbcVertexMap.h>

// Helix Hough includes
#ifndef __CINT__
#include <SimpleHit3D.h>
#include <SimpleTrack3D.h>
#include <VertexFinder.h> 
#include <sPHENIXTracker.h>
#endif

// ROOT includes
#include <TVector3.h>

// standard includes
#include <vector>
#include <map>
#include <float.h>

// g4hough includes
#include "SvtxTrackState.h"

// forward declarations
class PHCompositeNode;
class SvtxClusterMap;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class PHG4InEvent;
class PHG4HitContainer;

/// \class PHG4HoughTransform
///
/// \brief A fun4all implementation of Alan's Hough Transform
///
/// This module run Alan's Hough Transform and quick vertex finding
/// on the SvxClusterList of the event. It will produce both
/// SvxTrack and SvxSegments for the time being.
///
class PHG4HoughTransform : public SubsysReco {

public:
 
  PHG4HoughTransform(unsigned int seed_layers = 4,
		     unsigned int req_seed = 4,
		     const std::string &name = "PHG4HoughTransform");
  virtual ~PHG4HoughTransform() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  /// set verbosity
  void Verbosity(int verb) {verbosity = verb;}

  /// external handle for projecting tracks into the calorimetry
  static void projectToRadius(const SvtxTrack* track,
			      double magfield, // in Tesla
			      double radius,   // in cm
			      std::vector<double>& intersection);

  /// external handle for state objects
  static void projectToRadius(const SvtxTrackState* state,
			      int charge,
			      double magfield, // in Tesla
			      double radius,   // in cm
			      std::vector<double>& intersection);

  float get_mag_field() const          {return _magField;}
  void  set_mag_field(float magField) {_magField = magField;}
  
  /// set option to produce initial vertex for further tracking
  void set_use_vertex(bool b) {_use_vertex = b;}
  /// fetch option to produce initial vertex for further tracking
  bool get_use_vertex() {return _use_vertex;}

  /// set the tracking chi2 for initial vertex finding
  void set_chi2_cut_init(double chi2_cut) {_chi2_cut_init = chi2_cut;}
  /// get the tracking chi2 cut for initial vertex finding
  double get_chi2_cut_init() {return _chi2_cut_init;}

  /// set the tracking pt-dependent chi2 cut for fast fit, cut = min(par0 + par1 / pt, max)
  void set_chi2_cut_fast(double cut_par0,
			 double cut_par1 = 0.0,
			 double cut_max = FLT_MAX) {
    _chi2_cut_fast_par0 = cut_par0;
    _chi2_cut_fast_par1 = cut_par1;
    _chi2_cut_fast_max = cut_max;
  }

  /// set the tracking chi2 cut for full fit 
  void set_chi2_cut_full(double chi2_cut) {_chi2_cut_full = chi2_cut;}
  /// get the tracking chi2 cut for full fit
  double get_chi2_cut_full() {return _chi2_cut_full;}

  /// set early combination-land chi2 cut(?)
  void set_ca_chi2_cut(double chi2_cut) {_ca_chi2_cut = chi2_cut;}
  /// get early combination-land chi2 cut(?)
  double get_ca_chi2_cut() {return _ca_chi2_cut;}

  /// set early curvature cut between hits, lower values are more open
  void set_cos_angle_cut(double cos_angle_cut) {_cos_angle_cut = cos_angle_cut;}
  /// get early curvature cut between hits, lower values are more open
  double get_cos_angle_cut() {return _cos_angle_cut;}

  void setInitialResMultiplier(int beta) {_beta = beta;}
  void setFullResMultiplier(int lambda) {_lambda = lambda;}
  
  /// set the minimum pT to try to find during full tracking
  void set_min_pT(float PT) {_min_pT=PT;}
  /// set the minimum pT to try to find during initial vertex finding tracking
  void set_min_pT_init(float PT) {_min_pT_init=PT;}
  
  /// radiation length per layer, sequential layer indexes required here
  void set_material(int layer, float value);

  /// set internal ghost rejection
  void setRejectGhosts(bool rg) {_reject_ghosts = rg;}
  /// set for internal hit rejection
  void setRemoveHits(bool rh) {_remove_hits = rh;}

  /// use the cell size as cluster size instead of value stored on cluster
  void setUseCellSize(bool use_cell_size) {_use_cell_size = use_cell_size;}

  /// limit the maximum error reported by cluster (in number of cell units)
  void setMaxClusterError(float max_cluster_error) {_max_cluster_error = max_cluster_error;}

  /// adjusts the rate of zooming
  void setBinScale(float scale) {_bin_scale = scale;}
  /// adjusts the rate of zooming
  void setZBinScale(float scale) {_z_bin_scale = scale;}

  /// turn on DCA limitation
  void setCutOnDCA(bool cod) {_cut_on_dca = cod;}
  /// sets an upper limit on X-Y DCA
  void setDCACut(float dcut) {_dca_cut = dcut;}
  /// sets an upper limit on Z DCA
  void setDCAZCut(float dzcut) {_dcaz_cut = dzcut;}

  /// adjust the fit pt by a recalibration factor (constant B versus real mag
  /// field)
  void setPtRescaleFactor(float pt_rescale) {_pt_rescale = pt_rescale;}

  /// adjust the relative voting error scale w.r.t. the cluster size
  void setVoteErrorScale(unsigned int layer, float scale) {
    if (scale > 0.0) {
      _vote_error_scale.at(layer) = scale;
    } else {
      std::cout << "PHG4HoughTransform::setVoteErrorScale : scale must be "
                   "greater than zero ... doing nothing"
                << std::endl;
    }
  }
  /// adjust the relative fit error scale w.r.t. the cluster size
  void setFitErrorScale(unsigned int layer, float scale) {
    if (scale > 0.0) {
      _fit_error_scale.at(layer) = scale;
    } else {
      std::cout << "PHG4HoughTransform::setFitErrorScale : scale must be "
                   "greater than zero ... doing nothing"
                << std::endl;
    }
  }

#ifndef __CINT__
private:
  
  //--------------
  // InitRun Calls
  //--------------

  int GetNodes(PHCompositeNode *topNode);
  int CreateNodes(PHCompositeNode *topNode);
  int InitializeGeometry(PHCompositeNode *topNode);

  /// code to setup seed tracking objects
  int setup_seed_tracker_objects();
  
  /// code to setup initial vertexing tracker
  int setup_initial_tracker_object();
  
  /// code to setup full tracking object
  int setup_tracker_object();
  
  //--------------------
  // Process Event Calls
  //--------------------

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
  
  /// code to perform the final tracking and vertexing
  int full_tracking_and_vertexing();
  
  /// code to translate back to the SVTX universe
  int export_output();

  //------------------
  // Subfunction Calls
  //------------------

  /// convert from inverse curvature to momentum
  float kappaToPt(float kappa);
  
  /// convert from momentum to inverse curvature
  float ptToKappa(float pt);

  /// translate the clusters, tracks, and vertex from one origin to another
  void shift_coordinate_system(double dx, double dy, double dz);
  
  /// helper function for projection code
  static bool circle_line_intersections(double x0, double y0, double r0,
					double x1, double y1, double vx1, double vy1,
					std::set<std::vector<double> >* points);
  /// helper function for projection code
  static bool circle_circle_intersections(double x0, double y0, double r0,
					  double x1, double y1, double r1,
					  std::set<std::vector<double> >* points);

  BbcVertexMap* _bbc_vertexes;

  bool new_dca_nbin, new_z_z0, new_circle_dca, new_circle_kappa;
  
  bool _use_vertex;
  int _beta, _lambda; ///< resolution tuning parameters 

  double _chi2_cut_init;            ///< fit quality chisq/dof for initial track finding
  double _chi2_cut_fast_par0;       ///< fit quality chisq/dof for fast fit track fitting
  double _chi2_cut_fast_par1;       ///< fit quality chisq/dof for fast fit track fitting
  double _chi2_cut_fast_max;        ///< fit quality chisq/dof for fast fit track fitting
  double _chi2_cut_full;            ///< fit quality chisq/dof for kalman track fitting
  double _ca_chi2_cut;
  double _cos_angle_cut;
  unsigned int _max_hits_init;
  unsigned int _min_hits_init;
  unsigned int _maxtracks;
  unsigned int _max_hits;
  unsigned int _min_hits;

  int _nlayers;                        ///< number of detector layers
  std::vector<float> _radii;           ///< radial distance of each layer (cm)
  std::vector<float> _smear_xy_layer;  ///< detector hit resolution in phi (cm)
  std::vector<float> _smear_z_layer;   ///< detector hit resolution in z (cm)
  std::vector<float> _material;        ///< material at each layer in rad. lengths

  // object storage                                                                                                     
  std::vector<SimpleHit3D> _clusters;    ///< working array of clusters
  std::vector<SimpleTrack3D> _tracks;    ///< working array of tracks
  std::vector<double> _track_errors;     ///< working array of track chisq
  std::vector<Eigen::Matrix<float,5,5> > _track_covars; ///< working array of track covariances
  std::vector<float> _vertex;            ///< working array for collision vertex

  // track finding routines                                                                                             
  sPHENIXTracker *_tracker;    // finds full tracks
  sPHENIXTracker* _tracker_vertex; // finds a subset of tracks for initial vertex-finding
  sPHENIXTracker* _tracker_etap_seed; 
  sPHENIXTracker* _tracker_etam_seed;
  
  VertexFinder _vertexFinder; ///< vertex finding object

  float _magField; ///< in Tesla
  static float _cmToGeV;  ///< radius of curvature conversion (radius of curvature for a 1 GeV/c particle in 1 Tesla is 333.6 cm)

  std::vector<double> _radius_vec;
  
  SvtxClusterMap* _g4clusters;
  SvtxTrackMap* _g4tracks;
  SvtxVertexMap* _g4vertexes;
  PHG4InEvent* _inevent;

  float _min_pT;
  float _min_pT_init;
  
  PHG4HitContainer* _ghitlist;
  
  unsigned int _seed_layers, _req_seed;
  
  std::map<int, float> _user_material;
  
  bool _reject_ghosts;
  bool _remove_hits;
  bool _use_cell_size;

  float _max_cluster_error;
  
  float _bin_scale;
  float _z_bin_scale;
  
  bool _cut_on_dca;
  float _dca_cut;
  float _dcaz_cut;

  float _pt_rescale;
  
  std::vector<float> _fit_error_scale;
  std::vector<float> _vote_error_scale;

  /// recorded layer indexes to internal sequential indexes
  std::map<int,unsigned int> _layer_ilayer_map; 

#endif // __CINT__
};

#endif // __PHG4HOUGHTRANSFORM_H__
