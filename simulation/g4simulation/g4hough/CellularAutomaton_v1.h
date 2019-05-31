#ifndef G4HOUGH_CELLULARAUTOMATONV1_H
#define G4HOUGH_CELLULARAUTOMATONV1_H

#include "CellularAutomaton.h"

#include "Cluster3D.h"          // for Cluster3D

#if !defined(__CINT__) || defined(__CLING__)
#include "HelixTrackState.h"    // for HelixTrackState
#endif

#include <iostream>             // for cout, ostream
#include <map>                  // for map
#include <set>                  // for set
#include <vector>               // for vector

class HelixHoughSpace;
class HelixKalmanFilter;

class TrackSegment {
 public:
  TrackSegment()
      : chi2(0.), ux(0.), uy(0.), kappa(0.), dkappa(0.), seed(0), n_hits(0) {}
  ~TrackSegment() {}

  float chi2;
  float ux, uy;
  float kappa;
  float dkappa;
  int helicity;
  std::vector<unsigned int> hits; // cluster_ids
  unsigned int seed;
  unsigned int bin;
  unsigned int n_hits;
};



class CellularAutomaton_v1 : public CellularAutomaton {

 public:

	CellularAutomaton_v1(std::vector<Track3D>& input_tracks, std::vector<float>& radius, std::vector<float>& material);
//	CellularAutomaton_v1(const CellularAutomaton_v1& cellular_automaton);
	virtual ~CellularAutomaton_v1() {};

	// The "standard PHObject response" functions...
	void identify(std::ostream &os=std::cout) const {};
	void Reset();
	int  isValid() const {return 1;}
	CellularAutomaton* Clone() const {return new CellularAutomaton_v1(*this);}

	void set_hough_space(HelixHoughSpace* hough_space);
	void set_mag_field(float mag_field);
	void set_pt_rescale(float pt_rescale);
        void set_n_layers(unsigned int n_layers) {nlayers = n_layers;}
        void set_required_layers(unsigned int req_layers) {rlayers = req_layers;}
	void set_ca_chi2(float chi2_cut){ ca_chi2_cut = chi2_cut;}
	void set_ca_chi2_layer(float chi2_layer_cut){ca_chi2_layer_cut = chi2_layer_cut;}
	void set_ca_phi_cut(float phi_cut){ca_phi_cut = phi_cut;}
	void set_ca_z_cut(float z_cut) {ca_z_cut = z_cut;}
	void set_ca_dcaxy_cut(float dcaxy_cut){ca_dcaxy_cut = dcaxy_cut;}
	void set_propagate_forward(bool fwd) {forward = fwd;}
	void set_remove_hits(bool remove){remove_hits = remove;}
	void set_remove_inner_hits(bool remove_inner){remove_inner_hits = remove_inner;}
	void set_require_inner_hits(bool require_inner){require_inner_hits = require_inner;}
	void set_triplet_mode(bool mod) {triplet_mode = mod;}
	void set_seeding_mode(bool mod) {seeding_mode = mod;}
	void set_hits_map(std::map<unsigned int, Cluster3D>& hits_map){_hits_map = hits_map;}

	int run(std::vector<Track3D>& output_tracks, std::vector<HelixTrackState>& output_track_states, std::map<unsigned int, bool>& hits_used);	


 private:

	void set_detector_radii(std::vector<float>& radii);
	void set_detector_materials(std::vector<float>& materials);
	void set_input_tracks(std::vector<Track3D>& input_tracks);
	void set_cylinder_kalman();

	int init();
	int process_tracks();
	int process_single_track(Track3D& track);
	int process_single_triplet(Track3D& track);
	int get_ca_tracks(std::vector<Track3D>& output_tracks, std::vector<HelixTrackState>& output_track_states);	

	int calculate_kappa_tangents(                        
			float x1, float y1, float z1, float x2, float y2, float z2, 
                        float x3, float y3, float z3,
                        float dx1, float dy1, float dz1, float dx2, float dy2, float dz2,
                        float dx3, float dy3, float dz3, 
                        float& kappa, float& dkappa,
                        float& ux_mid, float& uy_mid, float& ux_end, float& uy_end,
                        float& dzdl_1, float& dzdl_2, float& ddzdl_1, float& ddzdl_2);
	int calculate_kappa_tangents(
                        float x1, float y1, float z1, float x2, float y2, float z2,
                        float x3, float y3, float z3,
                        float dx1, float dy1, float dz1, float dx2, float dy2, float dz2,
                        float dx3, float dy3, float dz3,
                        float& kappa, float& dkappa,
                        float& ux_mid, float& uy_mid, float& ux_end, float& uy_end,
                        float& dzdl_1, float& dzdl_2, float& ddzdl_1, float& ddzdl_2,
                        float ca_sin_ang_cut, float inv_cos_ang_diff,
                        float cur_kappa, float cur_dkappa, float cur_ux, float cur_uy, 
                        float cur_chi2, float& chi2);	

	float kappa_to_pt(float kappa);
	float shift_phi_range(float phi);

	HelixHoughSpace* _hough_space;
	HelixKalmanFilter* _kalman;
	std::vector<Track3D> in_tracks;
	std::vector<Track3D> ca_tracks;
	std::vector<HelixTrackState> ca_track_states;

	std::vector<unsigned int> temp_combo;
  	std::set<std::vector<unsigned int> > combos;
        std::vector<std::vector<Cluster3D> > layer_sorted;

	unsigned int nlayers; // number of layers for seeding
	unsigned int rlayers; // number of layers with hits required for seeding
	unsigned int allowed_missing_inner_hits;

        float ca_cos_ang_cut;
        float ca_chi2_cut;
	float ca_chi2_layer_cut;
	float ca_phi_cut;
	float ca_z_cut;
	float ca_dcaxy_cut;	

	float fast_chi2_cut_max;
	float fast_chi2_cut_par0;
	float fast_chi2_cut_par1;


	float _pt_rescale;
	float _mag_field;

	std::vector<float> _detector_radii;
	std::vector<float> _detector_scatter;
	std::vector<float> _detector_materials;
	std::vector<float> _integrated_scatter;

	std::map<unsigned int, bool> _hits_used;
	std::map<unsigned int, Cluster3D> _hits_map;

	double CAtime;
	double KALtime;	
	bool forward;
	bool remove_hits;
	bool remove_inner_hits;
	bool require_inner_hits;
	bool triplet_mode;
	bool seeding_mode;

	ClassDef(CellularAutomaton_v1,1)
};

#endif
