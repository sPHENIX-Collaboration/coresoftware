/*!
 *  \file PHG4KalmanPatRec.h
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

#ifndef __PHG4KALMANPATREC_H__
#define __PHG4KALMANPATREC_H__

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHTimeServer.h>
#include <phool/PHTimer.h>
#include <g4bbc/BbcVertexMap.h>

// Helix Hough includes
#ifndef __CINT__
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/VertexFinder.h> 
#include <HelixHough/sPHENIXSeedFinder.h>
#endif

// standard includes
#include <vector>
#include <map>
#include <memory>
#include <float.h>

// PHGenFit
#include <phgenfit/Fitter.h>

// g4hough includes
#include "SvtxTrackState.h"

// forward declarations
class PHCompositeNode;
class SvtxClusterMap;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class SvtxVertex;
class PHG4HitContainer;

namespace PHGenFit {
class Fitter;
class Track;
} /* namespace PHGenFit */

namespace genfit {
class GFRaveVertexFactory;
} /* namespace genfit */

/// \class PHG4KalmanPatRec
///
/// \brief A fun4all implementation of Alan's Hough Transform
///
/// This module run Alan's Hough Transform and quick vertex finding
/// on the SvxClusterList of the event. It will produce both
/// SvxTrack and SvxSegments for the time being.
///
class PHG4KalmanPatRec: public SubsysReco {

public:

	PHG4KalmanPatRec(unsigned int nlayers = 7, unsigned int min_nlayers = 7,
			const std::string &name = "PHG4KalmanPatRec");
	virtual ~PHG4KalmanPatRec() {
	}

	int Init(PHCompositeNode *topNode);
	int InitRun(PHCompositeNode *topNode);
	int process_event(PHCompositeNode *topNode);
	int End(PHCompositeNode *topNode);

	/// external handle for projecting tracks into the calorimetry
	static void projectToRadius(const SvtxTrack* track, double magfield, // in Tesla
			double radius,   // in cm
			std::vector<double>& intersection);

	/// external handle for state objects
	static void projectToRadius(const SvtxTrackState* state, int charge,
			double magfield, // in Tesla
			double radius,   // in cm
			std::vector<double>& intersection);

	void set_mag_field(float magField) {
		_magField = magField;
	}
	float get_mag_field() const {
		return _magField;
	}

	/// set the tracking pt-dependent chi2 cut for fast fit, cut = min(par0 + par1 / pt, max)
	void set_chi2_cut_fast(double cut_par0, double cut_par1 = 0.0,
			double cut_max = FLT_MAX) {
		_chi2_cut_fast_par0 = cut_par0;
		_chi2_cut_fast_par1 = cut_par1;
		_chi2_cut_fast_max = cut_max;
	}

	/// set the tracking chi2 cut for full fit
	void set_chi2_cut_full(double chi2_cut) {
		_chi2_cut_full = chi2_cut;
	}
	/// get the tracking chi2 cut for full fit
	double get_chi2_cut_full() const {
		return _chi2_cut_full;
	}

	/// set early combination-land chi2 cut(?)
	void set_ca_chi2_cut(double chi2_cut) {
		_ca_chi2_cut = chi2_cut;
	}
	/// get early combination-land chi2 cut(?)
	double get_ca_chi2_cut() const {
		return _ca_chi2_cut;
	}

	/// set early curvature cut between hits, lower values are more open
	void set_cos_angle_cut(double cos_angle_cut) {
		_cos_angle_cut = cos_angle_cut;
	}
	/// get early curvature cut between hits, lower values are more open
	double get_cos_angle_cut() const {
		return _cos_angle_cut;
	}

	void set_min_pT(float pt) {
		_min_pt = pt;
	}

	/// set the z0 search window
	void set_z0_range(float min_z0, float max_z0) {
		_min_z0 = min_z0;
		_max_z0 = max_z0;
	}

	/// set the r search window
	void set_r_max(float max_r) {
		_max_r = max_r;
	}

	/// radiation length per layer, sequential layer indexes required here
	void set_material(int layer, float value);

	/// set internal ghost rejection
	void setRejectGhosts(bool rg) {
		_reject_ghosts = rg;
	}
	/// set for internal hit rejection
	void setRemoveHits(bool rh) {
		_remove_hits = rh;
	}

	/// adjusts the rate of zooming
	void setBinScale(float scale) {
		_bin_scale = scale;
	}
	/// adjusts the rate of zooming
	void setZBinScale(float scale) {
		_z_bin_scale = scale;
	}

	/// turn on DCA limitation
	void setCutOnDCA(bool cod) {
		_cut_on_dca = cod;
	}
	/// sets an upper limit on X-Y DCA
	void setDCACut(float dcut) {
		_dcaxy_cut = dcut;
	}
	/// sets an upper limit on Z DCA
	void setDCAZCut(float dzcut) {
		_dcaz_cut = dzcut;
	}

	/// adjust the fit pt by a recalibration factor (constant B versus real mag
	/// field)
	void setPtRescaleFactor(float pt_rescale) {
		_pt_rescale = pt_rescale;
	}

	/// adjust the relative voting error scale w.r.t. the cluster size
	void setVoteErrorScale(unsigned int layer, float scale) {
		if (scale > 0.0) {
			_vote_error_scale.at(layer) = scale;
		} else {
			std::cout << "PHG4KalmanPatRec::setVoteErrorScale : scale must be "
					"greater than zero ... doing nothing" << std::endl;
		}
	}
	/// adjust the relative fit error scale w.r.t. the cluster size
	void setFitErrorScale(unsigned int layer, float scale) {
		if (scale > 0.0) {
			_fit_error_scale.at(layer) = scale;
		} else {
			std::cout << "PHG4KalmanPatRec::setFitErrorScale : scale must be "
					"greater than zero ... doing nothing" << std::endl;
		}
	}

	//---deprecated---------------------------------------------------------------

	/// set option to produce initial vertex for further tracking
	void set_use_vertex(bool b) {
	}

	void setInitialResMultiplier(int beta) {
	}
	void setFullResMultiplier(int lambda) {
	}

	/// set the minimum pT to try to find during initial vertex finding tracking
	void set_min_pT_init(float PT) {
	}

	/// limit the maximum error reported by cluster (in number of cell units)
	void setMaxClusterError(float max_cluster_error) {
	}

	/// use the cell size as cluster size instead of value stored on cluster
	void setUseCellSize(bool use_cell_size) {
	}

	/// set the tracking chi2 for initial vertex finding
	void set_chi2_cut_init(double chi2_cut) {
	}

	const std::vector<int>& get_seeding_layer() const {
		return _seeding_layer;
	}

//	void set_seeding_layer(const std::vector<int>& seedingLayer) {
//		_seeding_layer = seedingLayer;
//	}

	void set_seeding_layer(const int* seedingLayer, const int n) {
		_seeding_layer.assign(seedingLayer, seedingLayer + n);
	}

	float get_search_win_rphi() const {
		return _search_win_rphi;
	}

	void set_search_win_rphi(float searchWinMultiplier) {
		_search_win_rphi = searchWinMultiplier;
	}


	bool is_do_evt_display() const {
		return _do_evt_display;
	}

	void set_do_evt_display(bool doEvtDisplay) {
		_do_evt_display = doEvtDisplay;
	}

	const std::string& get_track_fitting_alg_name() const {
		return _track_fitting_alg_name;
	}

	void set_track_fitting_alg_name(const std::string& trackFittingAlgName) {
		_track_fitting_alg_name = trackFittingAlgName;
	}

	bool is_seeding_only_mode() const {
		return _seeding_only_mode;
	}

	void set_seeding_only_mode(bool seedingOnlyMode) {
		_seeding_only_mode = seedingOnlyMode;
	}

	float get_max_merging_deta() const {
		return _max_merging_deta;
	}

	void set_max_merging_deta(float maxMergingDeta) {
		_max_merging_deta = maxMergingDeta;
	}

	float get_max_merging_dphi() const {
		return _max_merging_dphi;
	}

	void set_max_merging_dphi(float maxMergingDphi) {
		_max_merging_dphi = maxMergingDphi;
	}

	float get_max_merging_dr() const {
		return _max_merging_dr;
	}

	void set_max_merging_dr(float maxMergingDr) {
		_max_merging_dr = maxMergingDr;
	}

	float get_max_merging_dz() const {
		return _max_merging_dz;
	}

	void set_max_merging_dz(float maxMergingDz) {
		_max_merging_dz = maxMergingDz;
	}


	float get_search_win_z() const {
		return _search_win_z;
	}

	void set_search_win_z(float searchWinZ) {
		_search_win_z = searchWinZ;
	}

	float get_max_incr_chi2() const {
		return _max_incr_chi2;
	}

	void set_max_incr_chi2(float maxIncrChi2) {
		_max_incr_chi2 = maxIncrChi2;
	}

	unsigned int get_max_consecutive_missing_layer() const {
		return _max_consecutive_missing_layer;
	}

	void set_max_consecutive_missing_layer(
			unsigned int maxConsecutiveMissingLayer) {
		_max_consecutive_missing_layer = maxConsecutiveMissingLayer;
	}

	unsigned int get_min_good_track_hits() const {
		return _min_good_track_hits;
	}

	void set_min_good_track_hits(unsigned int minGoodTrackHits) {
		_min_good_track_hits = minGoodTrackHits;
	}

	unsigned int get_max_share_hits() const {
		return _max_share_hits;
	}

	void set_max_share_hits(unsigned int maxShareHits) {
		_max_share_hits = maxShareHits;
	}

	float get_max_splitting_chi2() const {
		return _max_splitting_chi2;
	}

	void set_max_splitting_chi2(float maxSplittingChi2) {
		_max_splitting_chi2 = maxSplittingChi2;
	}

#ifndef __CINT__

private:

	//--------------
	// InitRun Calls
	//--------------

	/// create new node output pointers
	int CreateNodes(PHCompositeNode *topNode);

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
	int GetNodes(PHCompositeNode *topNode);

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
	int full_track_seeding();

	/// code to translate back to the SVTX universe
	int export_output();

	//!
	int CleanupSeeds();

	//!
	int FullTrackFitting(PHCompositeNode* topNode);

	//!
	int ExportOutput();

	//--------------------
	//
	//--------------------

	//! FullTrackFitting Call.
	int BuildLayerZPhiHitMap();

	//! layer: 7 bits, z: 11 bits, phi: 14 bits
	unsigned int encode_cluster_index(const unsigned int layer, const float z, const float rphi);

	unsigned int encode_cluster_index(const unsigned int layer, const unsigned int iz, const unsigned int irphi);

	//! FullTrackFitting Call.
	int SimpleTrack3DToPHGenFitTracks(PHCompositeNode* topNode, unsigned int itrack);
	int TrackPropPatRec(PHCompositeNode* topNode, const int iPHGenFitTrack, std::shared_ptr<PHGenFit::Track> track);

	//! TrackPropPatRec Call.
	std::vector<unsigned int> SearchHitsNearBy (const unsigned int layer, const float z_center, const float phi_center, const float z_window, const float phi_window);


	//! ExportOutput Call. Make SvtxTrack from PHGenFit::Track and set of clusters
	//std::shared_ptr<SvtxTrack> MakeSvtxTrack(const int genfit_track_ID, const SvtxVertex * vertex = NULL);
	int OutputPHGenFitTrack(PHCompositeNode* topNode, std::map<int, std::shared_ptr<PHGenFit::Track>>::iterator);

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

	/// helper function for projection code
	static bool circle_line_intersections(double x0, double y0, double r0,
			double x1, double y1, double vx1, double vy1,
			std::set<std::vector<double> >* points);

	/// helper function for projection code
	static bool circle_circle_intersections(double x0, double y0, double r0,
			double x1, double y1, double r1,
			std::set<std::vector<double> >* points);

	int _event;
	PHTimer *_t_seeding;
	PHTimer *_t_seeds_cleanup;
	PHTimer *_t_translate_to_PHGenFitTrack;
	PHTimer *_t_kalman_pat_rec;
	PHTimer *_t_search_clusters;
	PHTimer *_t_search_clusters_encoding;
	PHTimer *_t_search_clusters_map_iter;
	PHTimer *_t_track_propagation;
	PHTimer *_t_full_fitting;
	PHTimer *_t_output_io;

	std::vector<int> _seeding_layer; //layer numbers that are used for seeding

	unsigned int _nlayers;               ///< number of detector layers
	unsigned int _min_nlayers;     ///< minimum number of layers to make a track
	std::vector<float> _radii;           ///< radial distance of each layer (cm)
	std::vector<float> _material;    ///< material at each layer in rad. lengths
	std::map<int, float> _user_material; ///< material in user ladder indexes

	float _magField; ///< in Tesla

	bool _reject_ghosts;
	bool _remove_hits;

	float _min_pt;
	float _min_z0;
	float _max_z0;
	float _max_r;

	bool _cut_on_dca;
	float _dcaxy_cut;
	float _dcaz_cut;

	double _chi2_cut_fast_par0; ///< fit quality chisq/dof for fast fit track fitting
	double _chi2_cut_fast_par1; ///< fit quality chisq/dof for fast fit track fitting
	double _chi2_cut_fast_max; ///< fit quality chisq/dof for fast fit track fitting
	double _chi2_cut_full;   ///< fit quality chisq/dof for kalman track fitting
	double _ca_chi2_cut;              ///< initial combination cut?
	double _cos_angle_cut;          ///< curvature restriction on cluster combos

	float _bin_scale;
	float _z_bin_scale;

	unsigned int _min_combo_hits; ///< minimum hits to enter combination gun
	unsigned int _max_combo_hits; ///< maximum hits to enter combination gun

	float _pt_rescale;
	std::vector<float> _fit_error_scale;
	std::vector<float> _vote_error_scale;

	/// recorded layer indexes to internal sequential indexes
	std::map<int, unsigned int> _layer_ilayer_map;

	// object storage
	std::vector<SimpleHit3D> _clusters;    ///< working array of clusters
	std::vector<SimpleTrack3D> _tracks;    ///< working array of tracks
	std::vector<double> _track_errors;     ///< working array of track chisq
	std::vector<Eigen::Matrix<float, 5, 5> > _track_covars; ///< working array of track covariances
	std::vector<float> _vertex;          ///< working array for collision vertex

	// track finding routines
	sPHENIXSeedFinder *_tracker;           ///< finds full tracks
	sPHENIXSeedFinder* _tracker_vertex; ///< finds a subset of tracks for initial vertex-finding
	sPHENIXSeedFinder* _tracker_etap_seed; ///< finds a subset of tracks for the vertex guess
	sPHENIXSeedFinder* _tracker_etam_seed; ///< finds a subset of tracks for the vertex guess
	VertexFinder _vertexFinder;      ///< vertex finding object

	// node pointers
	BbcVertexMap* _bbc_vertexes;
	SvtxClusterMap* _g4clusters;
	SvtxTrackMap* _g4tracks;
	SvtxVertexMap* _g4vertexes;


	bool _seeding_only_mode;

	//! Cleanup Seeds
	float _max_merging_dphi;
	float _max_merging_deta;
	float _max_merging_dr;
	float _max_merging_dz;

	//! if two seeds have common hits more than this number, merge them
	unsigned int _max_share_hits;

	//GenFit Related Part
	//!
	std::string _mag_field_file_name;

	//! rescale mag field, modify the original mag field read in
	float _mag_field_re_scaling_factor;

	//! Switch to reverse Magnetic field
	bool _reverse_mag_field;

	PHGenFit::Fitter* _fitter;

	//! KalmanFitterRefTrack, KalmanFitter, DafSimple, DafRef
	//PHGenFit::Fitter::FitterType _track_fitting_alg_name;
	std::string _track_fitting_alg_name;

	int _primary_pid_guess;
	double _cut_min_pT;

	bool _do_evt_display;

	int _nlayers_all;
	std::map<int, unsigned int> _layer_ilayer_map_all;
	std::vector<float> _radii_all;
	float _search_win_rphi;
	float _search_win_z;
	std::map<int, float> _search_wins_rphi;
	std::map<int, float> _search_wins_z;

	//std::map<unsigned int, std::map<int, std::multimap<int, unsigned int>>> _layer_zID_phiID_cluserID;
	std::multimap<unsigned int,  unsigned int> _layer_zID_phiID_cluserID;

	float _half_max_z;
	float _half_max_rphi;
	float _layer_zID_phiID_cluserID_phiSize;
	float _layer_zID_phiID_cluserID_zSize;


	std::map<int, std::shared_ptr<PHGenFit::Track>> _trackID_PHGenFitTrack;
	//std::map<int, std::vector<unsigned int>> _trackID_clusterID;

	unsigned int _max_consecutive_missing_layer;
	float _max_incr_chi2;
	std::map<int, float> _max_incr_chi2s;

	float _max_splitting_chi2;

	unsigned int _min_good_track_hits;

#endif // __CINT__
};

#endif // __PHG4KALMANPATREC_H__
