#ifndef TRACKRECO_PHINITZVERTEXING_H
#define TRACKRECO_PHINITZVERTEXING_H

#include "PHInitVertexing.h"   // base class
#include "VertexFitter.h"


// Helix Hough includes
#if !defined(__CINT__) || defined(__CLING__)
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/HelixKalmanState.h>
#include <Eigen/Core>            // for Matrix
#endif

#include <trackbase/TrkrDefs.h>

// standard includes
#include <algorithm>             // for copy
#include <map>
#include <string>                // for string
#include <vector>


// forward declarations
class BbcVertexMap;

//class CellularAutomaton;
class HelixHoughBin;
class HelixHoughSpace;   
class HelixHoughFuncs;

class PHCompositeNode;
class PHTimer;

class TrkrClusterContainer;
class SvtxTrackMap;
class SvtxVertexMap;

class TFile;
class TH2D;
class TH3D;
class TNtuple;

///
class PHInitZVertexing: public PHInitVertexing {

public:

	PHInitZVertexing(unsigned int nlayers = 7, unsigned int min_nlayers = 7,
			const std::string &name = "PHInitZVertexing");
	virtual ~PHInitZVertexing() {
	}

	void set_file_name(const std::string &fname){_fname = fname;}
	void set_mag_field(float mag_field) {
		_mag_field = mag_field;
	}
	float get_mag_field() const {
		return _mag_field;
	}

	void set_min_pT(float pt) {
		_min_pt = pt;
	}

	void set_max_kappa(float kappa_max);
	void set_min_d(float d_min) {_min_d = d_min;}
	void set_max_d(float d_max) {_max_d = d_max;}
	void set_min_dzdl(float dzdl_min) {_min_dzdl = dzdl_min; }
	void set_max_dzdl(float dzdl_max) {_max_dzdl = dzdl_max;}
	void set_min_z0(float z0_min) { _min_z0 = z0_min;}
	void set_max_z0(float z0_max) { _max_z0 = z0_max;}
	float get_max_kappa() const {return _max_kappa;}
	float get_min_z0() const {return _min_z0;}
	float get_max_z0() const {return _max_z0;}
	
	/// radiation length per layer, sequential layer indexes required here
	void set_material(int layer, float value);

	/// turn on DCA limitation
	void setCutOnDCA(bool cod) {
		_cut_on_dca = cod;
	}
	/// sets an upper limit on X-Y DCA
	void set_dcaxy_cut(float dcut) {
		_dcaxy_cut = dcut;
	}
	/// sets an upper limit on Z DCA
	void set_dcaz_cut(float dzcut) {
		_dcaz_cut = dzcut;
	}
	/// sets a chi2 cut 
	void set_chi2_cut(float chi2_cut) {
		_ca_chi2 = chi2_cut;
	}
	
	void set_ca_z_cut(float z_cut) {
		_ca_z_cut = z_cut;
	}


	/// adjust the fit pt by a recalibration factor (constant B versus real mag
	/// field)
	void setPtRescaleFactor(float pt_rescale) {
		_pt_rescale = pt_rescale;
	}

	void set_ca_phi_cut(float phi_cut)
	{	
		_ca_phi_cut = phi_cut;
	}

	void set_mult_onebin(float mult)
	{
		_mult_onebin = mult;
	}

	void set_mult_twobins(float mult)
	{
		_mult_twobins = mult;
	}

	void set_mult_threebins(float mult)
	{
		_mult_threebins= mult;
	}

	void set_min_zvtx_tracks(unsigned int min_zvtx_tracks)
	{
		_min_zvtx_tracks = min_zvtx_tracks;
	}

	const std::vector<int>& get_seeding_layer() const {
		return _seeding_layer;
	}

	void set_seeding_layer(const int* seedingLayer, const int n) {
		_seeding_layer.assign(seedingLayer, seedingLayer + n);
	}


	void add_zoom(unsigned int n_kappa, unsigned int n_phi, unsigned int n_d, unsigned int n_dzdl, unsigned int n_z0);
	void set_nzooms() {nzooms = zooms_vec.size();}
	void reset_zooms() {zooms_vec.clear();}

#if !defined(__CINT__) || defined(__CLING__)

 protected:

  int Setup(PHCompositeNode *topNode);
  int Process(PHCompositeNode *topNode);
  int End();

private:

  int Init(PHCompositeNode *topNode);

	//--------------
	// InitRun Calls
	//--------------
	int create_nodes(PHCompositeNode *topNode);
	int initialize_geometry(PHCompositeNode *topNode);

	//--------------------
	// Process Event Calls
	//--------------------
	int get_nodes(PHCompositeNode *topNode);
	int translate_input(PHCompositeNode *topNode);

        void set_nbins(unsigned int izoom);
        void initialize_houghbin();
        void vote_z_init(unsigned int zoomlevel);

        void find_track_candidates_z_init(unsigned int zoomlevel);
        void vote_z(unsigned int zoomlevel);
        void find_track_candidates_z(unsigned int zoomlevel);
        void vote_xy(unsigned int zoomlevel);
        void find_track_candidates_xy(unsigned int zoomlevel);

        void prune_z(unsigned int zoomlevel);
        void prune_xy(unsigned int zoomlevel);

	void reset_hits_used();

	int export_output();


	//------------------
	// Subfunction Calls
	//------------------
	/// convert from inverse curvature to momentum
	float kappa_to_pt(float kappa);

	/// convert from momentum to inverse curvature
	float pt_to_kappa(float pt);

	void bins_to_SimpleTrack3D(std::vector<SimpleTrack3D>& temp_tracks, int imap, unsigned int zoomlevel);
	int build_triplets_to_SimpleTrack3D(std::vector<SimpleTrack3D>& new_tracks, bool forward);
	int cellular_automaton_zvtx_init(std::vector<SimpleTrack3D>& temp_tracks);
//	int cellular_automaton_zvtx_second(std::vector<SimpleTrack3D>& temp_tracks);
//	int cellular_automaton_zvtx_third(std::vector<SimpleTrack3D>& temp_tracks);
	int fit_vertex();	

	/// convert the covariance from HelixHough coords to x,y,z,px,py,pz
	void convertHelixCovarianceToEuclideanCovariance(float B, float phi,
			float d, float kappa, float z0, float dzdl,
			Eigen::Matrix<float, 5, 5> const& input,
			Eigen::Matrix<float, 6, 6>& output);

	void shift_coordinate_system(double dx, double dy, double dz);
	float shift_phi_range(float phi);

	int _event;
	PHTimer *_t_output_io;

	std::vector<int> _seeding_layer; //layer numbers that are used for seeding

	unsigned int _nlayers;               ///< number of detector layers
	unsigned int _min_nlayers;     ///< minimum number of layers to make a track
	unsigned int _ca_nlayers;
	
	float _ca_chi2;
	std::vector<float> _radii;           ///< radial distance of each layer (cm)
	std::vector<float> _material;    ///< material at each layer in rad. lengths
	std::map<int, float> _user_material; ///< material in user ladder indexes

	float _mag_field; ///< in Tesla

	bool _use_max_kappa;
	bool fill_multi_zvtx;

	float _min_pt;
	float _max_kappa;
        float _min_d;
        float _max_d;
	float _min_dzdl;
	float _max_dzdl;
        float _min_z0;
        float _max_z0;

	bool _cut_on_dca;
	float _dcaxy_cut;
	float _dcaz_cut;
	float _pt_rescale;
	float _ca_phi_cut;
	float _ca_z_cut;
	float _z_cut;
	float _mult_onebin;
	float _mult_twobins;
	float _mult_threebins;
	unsigned int _min_zvtx_tracks;


        unsigned int bin;

	unsigned int ik;
	unsigned int ip;
	unsigned int id;
	unsigned int il;
	unsigned int iz;

	unsigned int nkappa;
	unsigned int nphi;
	unsigned int nd;
	unsigned int ndzdl;
	unsigned int nz0;

	//unsigned int cluster_id;
	TrkrDefs::cluskey cluster_key;

	/// recorded layer indexes to internal sequential indexes
	std::map<int, unsigned int> _layer_ilayer_map;

	// object storage
	std::vector<SimpleTrack3D> _temp_tracks;
        std::vector<SimpleTrack3D> _tracks;    ///< working array of tracks
        std::vector<HelixKalmanState> _track_states;
	std::vector<double> _track_errors;     ///< working array of track chisq
	std::vector<Eigen::Matrix<float, 5, 5> > _track_covars; ///< working array of track covariances
	std::vector<float> _vertex;          ///< working array for collision vertex
	std::vector<std::vector<float> > _vertex_list;

	// node pointers
	BbcVertexMap* _bbc_vertexes;
	SvtxTrackMap* _trackmap;
	VertexFitter _vertex_finder;

  	HelixHoughSpace* _hough_space;
	HelixHoughFuncs* _hough_funcs;

//	CellularAutomaton* ca;

	TNtuple* _ntp_zvtx_by_event;
	TNtuple* _ntp_zvtx_by_track;
	TH2D* _z0_dzdl;
	TH2D* _kappa_phi;
	TH2D* _d_phi;
	TH3D* _kappa_d_phi;
	TFile* _ofile;
	TFile* _ofile2;

	std::string _fname;
	int _nlayers_all;
	std::map<int, unsigned int> _layer_ilayer_map_all;
	std::vector<float> _radii_all;
	std::vector<std::vector<unsigned int> > zooms_vec;

	std::map<unsigned int, bool> hits_used;
	std::map<unsigned int, SimpleHit3D>   hits_map;
	std::map<unsigned int, HelixHoughBin* > bins_map;
        std::map<unsigned int, HelixHoughBin* > bins_map_prev;
        std::map<unsigned int, HelixHoughBin* > bins_map_cur;
        std::map<unsigned int, HelixHoughBin* > bins_map_sel;

	std::map<unsigned int, unsigned int> kappa_phi_map;
	std::map<unsigned int, unsigned int> d_phi_map;
	std::map<unsigned int, unsigned int> kappa_d_phi_map;

	float hit[50][3];
	float hitpos3d[3];
	unsigned int nzooms;

	bool separate_helicity;
	int helicity;
	int n_vtx_tracks;
	
#endif // !defined(__CINT__) || defined(__CLING__)
};

#endif // TRACKRECO_PHG4INITZVERTEXING_H
