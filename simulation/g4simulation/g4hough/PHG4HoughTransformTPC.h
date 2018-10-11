#ifndef PHG4HOUGHTRANSFORMTPC_H__
#define PHG4HOUGHTRANSFORMTPC_H__

//===========================================================
/// \file PHG4HoughTransformTPC.h
/// \brief A fun4all implementation of Alan's Hough Transform
/// \author Matt Wysocki (copied from SvxHoughTransform)
/// go to https://www.phenix.bnl.gov/WWW/offline/wikioffline/index.php/SvxHoughTransform
/// \edited by Theo Koblesky to conform to new changes in HelixHough (1/10/2012)

// edited week of 11/14/2012 by Alan Dion to use new HelixHough
// classes

//===========================================================

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
//#include <fun4all/Fun4AllReturnCodes.h>
//#include <g4bbc/BbcVertexMap.h>

// Helix Hough includes
#ifndef __CINT__
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/VertexFinder.h> 
#endif

// standard includes
#include <float.h>
#include <map>
#include <set>
#include <vector>

// forward declarations
class BbcVertexMap;
class PHCompositeNode;
class sPHENIXTrackerTPC;
class SvtxClusterMap;
class SvtxTrack;
class SvtxTrackMap;
class SvtxTrackState;
class SvtxVertexMap;
class SimpleRecoEvent;
class TTree;

/// \class PHG4HoughTransformTPC
///
/// \brief A fun4all implementation of Alan's Hough Transform
///
/// This module run Alan's Hough Transform and quick vertex finding
/// on the SvxClusterList of the event. It will produce both
/// SvxTrack and SvxSegments for the time being.
///
class PHG4HoughTransformTPC : public SubsysReco {

public:
 
  PHG4HoughTransformTPC(unsigned int seed_layers = 4,
		     unsigned int req_seed = 4,
		     const std::string &name = "PHG4HoughTransformTPC");
  virtual ~PHG4HoughTransformTPC() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  /// external handle for projecting tracks into the calorimetry
  static void projectToRadius(const SvtxTrack* track,
			      double magfield, // in Tesla
			      double radius,   // in cm
			      std::vector<double>& intersection);

  static void projectToRadius(const SvtxTrackState* state,
			      int charge,
			      double magfield, // in Tesla
			      double radius,   // in cm
			      std::vector<double>& intersection);

  float get_mag_field() const          {return _magField;}
  void  set_mag_field(float magField) {_magField = magField;}

  /// set the tracking pt-dependent chi2 cut for fast fit, cut = min(par0 + par1 / pt, max)
  void set_chi2_cut_fast(double cut_par0,
			 double cut_par1 = 0.0,
			 double cut_max = FLT_MAX) {
    _chi2_cut_fast_par0 = cut_par0;
    _chi2_cut_fast_par1 = cut_par1;
    _chi2_cut_fast_max = cut_max;
  }

  /// set the initial tracking chi2 cut
  void set_chi2_cut_init(double chi2_cut_init) {_chi2_cut_init = chi2_cut_init;}
  /// get the initial tracking chi2 cut
  double get_chi2_cut_init() { return _chi2_cut_init;}
  /// set the tracking chi2 cut for full fit 
  void set_chi2_cut_full(double chi2_cut) { _chi2_cut_full = chi2_cut;}
  /// get the tracking chi2 cut for full fit
  double get_chi2_cut_full() { return _chi2_cut_full;}

  /// set early combination-land chi2 cut(?)
  void set_ca_chi2_cut(double chi2_cut) { _ca_chi2_cut = chi2_cut;}
  /// get early combination-land chi2 cut(?)
  double get_ca_chi2_cut() { return _ca_chi2_cut;}

  /// set early curvature cut between hits, lower values are more open
  void set_cos_angle_cut(double cos_angle_cut) { _cos_angle_cut = cos_angle_cut;}
  /// get early curvature cut between hits, lower values are more open
  double get_cos_angle_cut() { return _cos_angle_cut;}
  
  /// set the minimum pt to try to find during full tracking
  void set_min_pT(float pt){_min_pt=pt;}

  /// set the z0 search window 
  void set_z0_range(float min_z0, float max_z0) {
    _min_z0 = min_z0; _max_z0 = max_z0;
  }

  /// set the d search window
  void set_d_max(float max_d) {_max_d = max_d;}

  /// radiation length per layer, sequential layer indexes required here
  void set_material(int layer, float value);

  /// set use bbc
  void set_use_vertex(bool ub){_use_bbc = ub;}  

  /// set internal ghost rejection
  void setRejectGhosts(bool rg){_reject_ghosts = rg;}
  /// set for internal hit rejection
  void setRemoveHits(bool rh){_remove_hits = rh;}

  /// use the cell size as cluster size instead of value stored on cluster
  void setUseCellSize(bool use_cell_size) {_use_cell_size = use_cell_size;}

  /// limit the maximum error reported by cluster (in number of cell units)
  void setMaxClusterError(float max_cluster_error) {_max_cluster_error = max_cluster_error;} 

  /// adjusts the rate of zooming
  void setBinScale(float scale){_bin_scale = scale;}
  /// adjusts the rate of zooming
  void setZBinScale(float scale){_z_bin_scale = scale;}

  /// turn on DCA limitation
  void setCutOnDCA(bool cod){_cut_on_dca = cod;}
  /// sets an upper limit on X-Y DCA
  void setDCACut(float dcut){_dcaxy_cut = dcut;}
  /// sets an upper limit on Z DCA
  void setDCAZCut(float dzcut){_dcaz_cut = dzcut;}

  /// adjust the fit pt by a recalibration factor (constant B versus real mag field)
  void setPtRescaleFactor(float pt_rescale) {_pt_rescale = pt_rescale;}

  /// adjust the relative voting error scale w.r.t. the cluster size
  void setVoteErrorScale(unsigned int layer, float scale) {
    if (scale > 0.0) {
      _vote_error_scale.at(layer) = scale;
    } else{
      std::cout << "PHG4HoughTransformTPC::setVoteErrorScale : scale must be "
                   "greater than zero ... doing nothing"
		<< std::endl;
    }
  }
  /// adjust the relative fit error scale w.r.t. the cluster size
  void setFitErrorScale(unsigned int layer, float scale) {
    if(scale > 0.0) {
      _fit_error_scale.at(layer) = scale;
    } else{
      std::cout << "PHG4HoughTransformTPC::setFitErrorScale : scale must be "
                   "greater than zero ... doing nothing"
                << std::endl;
    }
  }

  void setWriteRecoTree(bool flag){_write_reco_tree=flag;}

  //---deprecated---------------------------------------------------------------
/*
  /// set option to produce initial vertex for further tracking
  void set_use_vertex(bool b) {_use_vertex = b;}
  /// fetch option to produce initial vertex for further tracking
  bool get_use_vertex() {return _use_vertex;}

  /// set the tracking chi2 for initial vertex finding
  void set_chi2_cut_init(double chi2_cut) { _chi2_cut_init = chi2_cut;}
  /// get the tracking chi2 cut for initial vertex finding
  double get_chi2_cut_init() { return _chi2_cut_init;}
  
  void setInitialResMultiplier(int beta){ _beta = beta;}
  void setFullResMultiplier(int lambda){ _lambda = lambda;}
  /// set the minimum pT to try to find during initial vertex finding tracking
  void set_min_pT_init(float PT){_min_pT_init=PT;}

  double _chi2_cut_init; ///< fit quality chisq/dof for initial track finding
  int _beta, _lambda; ///< resolution tuning parameters 
  float _min_pT_init;
*/

#ifndef __CINT__
 private:

  //--------------
  // InitRun Calls
  //--------------

  /// create new node output pointers
  int CreateNodes(PHCompositeNode *topNode);

  /// scan tracker geometry objects
  int InitializeGeometry(PHCompositeNode *topNode);

  /// setup seed tracking objects
  int setup_seed_tracker_objects();

  /// setup initial vertexing tracker
  int setup_initial_tracker_object();

  /// code to setup full tracking object
  int setup_tracker_object();

  //--------------------
  // Process Event Calls
  //--------------------

  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);

  /// translate into the HelixHough universe
  int translate_input();

  /// combine seed tracking vertex with BBCZ if available
  int fast_composite_seed();
   
  /// seed vertex from bbc
  int fast_vertex_from_bbc();

  /// seed vertex from initial tracking using a broad search window
  int initial_zvertex_finding();
 
  /// perform the final tracking and vertexing
  int full_tracking_and_vertexing();

  /// translate back to the SVTX universe
  int export_output();

  //------------------
  // Subfunction Calls
  //------------------
 
  /// convert from inverse curvature to momentum
  float kappaToPt(float kappa);
  /// convert from momentum to inverse curvature
  float ptToKappa(float pt);

  /// convert the covariance from HelixHough coords to x,y,z,px,py,pz
  void convertHelixCovarianceToEuclideanCovariance(
      float B, float phi, float d, float kappa, float z0, float dzdl,
      Eigen::Matrix<float, 5, 5> const& input,
      Eigen::Matrix<float, 6, 6>& output);

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
  
  int _nlayers;              ///< number of detector layers                                                         
  unsigned int _seed_layers, _req_seed;
  std::vector<float> _radii;          ///< radial distance of each layer (cm)                                           
  std::vector<float> _smear_xy_layer; ///< detector hit resolution in phi (cm)                                          
  std::vector<float> _smear_z_layer;  ///< detector hit resolution in z (cm)                 
  std::vector<float> _material;  ///< material at each layer in rad. lengths
  std::map<int, float> _user_material;

  float _magField; ///< in Tesla
  static float _cmToGeV;  ///< radius of curvature conversion (radius of curvature for a 1 GeV/c particle in 1 Tesla is 333.6 cm)
  
  bool _use_bbc;
  bool _reject_ghosts;
  bool _remove_hits;
  bool _use_cell_size;
  float _max_cluster_error;

  float _min_pt;
  float _min_z0;
  float _max_z0;
  float _max_d;

  bool _cut_on_dca;
  float _dcaxy_cut;
  float _dcaz_cut;

  double _chi2_cut_fast_par0;       ///< fit quality chisq/dof for fast fit track fitting
  double _chi2_cut_fast_par1;       ///< fit quality chisq/dof for fast fit track fitting
  double _chi2_cut_fast_max;        ///< fit quality chisq/dof for fast fit track fitting
  double _chi2_cut_init;	    ///< fil quality chisq/dof for initial track fitting
  double _chi2_cut_full;            ///< fit quality chisq/dof for kalman track fitting
  double _ca_chi2_cut;
  double _cos_angle_cut;

  float _bin_scale;
  float _z_bin_scale;

  unsigned int _min_combo_hits;
  unsigned int _max_combo_hits;

  unsigned int _min_vtx_hits;
  unsigned int _max_vtx_hits;

  float _pt_rescale;
  std::vector<float> _fit_error_scale;
  std::vector<float> _vote_error_scale;

  /// recorded layer indexes to internal sequential indexes
  std::map<int,unsigned int> _layer_ilayer_map;

  // object storage                                                                                                     
  std::vector<SimpleHit3D> _clusters_init; ///< working array of clusters                                                    
  std::vector<SimpleHit3D> _clusters;
  std::vector<SimpleTrack3D> _tracks; ///< working array of tracks
  std::vector<double> _track_errors;     ///< working array of track chisq                                                      
  std::vector<Eigen::Matrix<float,5,5> > _track_covars; ///< working array of track covariances
  std::vector<float> _vertex;         ///< working array for collision vertex                                           

  // track finding routines                                                                                             
  sPHENIXTrackerTPC *_tracker;    // finds full tracks  
  sPHENIXTrackerTPC* _tracker_etap_seed; ///< finds a subset of tracks for the vertex guess
  sPHENIXTrackerTPC* _tracker_etam_seed; ///< finds a subset of tracks for the vertex guess
  VertexFinder _vertexFinder; ///< vertex finding object

  BbcVertexMap* _bbc_vertexes; 
  SvtxClusterMap* _g4clusters;
  SvtxTrackMap* _g4tracks;
  SvtxVertexMap* _g4vertexes;

  PHTimeServer::timer _timer;
  PHTimeServer::timer _timer_initial_hough;
  
  bool _write_reco_tree;
  TTree* _reco_tree;
  SimpleRecoEvent* _recoevent;

#endif // __CINT__
};

#endif // __SVXHOUGHTRANSFORMTPC_H__
