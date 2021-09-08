#ifndef G4HOUGH_CELLULARAUTOMATONV1_H
#define G4HOUGH_CELLULARAUTOMATONV1_H

#include "CellularAutomaton.h"

#include <HelixHough/SimpleHit3D.h>       
#include <HelixHough/SimpleTrack3D.h>     // for SimpleTrack3D

#include <HelixHough/HelixKalmanState.h>    // for HelixKalmanState

#include <iostream>             // for cout, ostream
#include <map>                  // for map
#include <set>                  // for set
#include <vector>               // for vector

class HelixHoughSpace;
class HelixKalmanFilter;

class TrackSegment {
 public:
  TrackSegment()
      : chi2(0.), ux(0.), uy(0.), kappa(0.), dkappa(0.), helicity(0), seed(0), bin(0), n_hits(0) {}
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

	CellularAutomaton_v1(std::vector<SimpleTrack3D>& input_tracks, std::vector<float>& radius, std::vector<float>& material);
//	CellularAutomaton_v1(const CellularAutomaton_v1& cellular_automaton);
	~CellularAutomaton_v1() override {};

	// The "standard PHObject response" functions...
	void identify(std::ostream &os=std::cout) const override {};
	void Reset() override;
	int  isValid() const override {return 1;}
	CellularAutomaton* Clone() const override {return new CellularAutomaton_v1(*this);}

	void set_hough_space(HelixHoughSpace* hough_space) override;
	void set_mag_field(float mag_field) override;
	void set_pt_rescale(float pt_rescale) override;
        void set_n_layers(unsigned int n_layers) override {nlayers = n_layers;}
        void set_required_layers(unsigned int req_layers) override {rlayers = req_layers;}
	void set_ca_chi2(float chi2_cut) override{ ca_chi2_cut = chi2_cut;}
	void set_ca_chi2_layer(float chi2_layer_cut) override{ca_chi2_layer_cut = chi2_layer_cut;}
	void set_ca_phi_cut(float phi_cut) override{ca_phi_cut = phi_cut;}
	void set_ca_z_cut(float z_cut) override {ca_z_cut = z_cut;}
	void set_ca_dcaxy_cut(float dcaxy_cut) override{ca_dcaxy_cut = dcaxy_cut;}
	void set_propagate_forward(bool fwd) override {forward = fwd;}
	void set_remove_hits(bool remove) override{remove_hits = remove;}
	void set_remove_inner_hits(bool remove_inner) override{remove_inner_hits = remove_inner;}
	void set_require_inner_hits(bool require_inner) override{require_inner_hits = require_inner;}
	void set_triplet_mode(bool mod) override {triplet_mode = mod;}
	void set_seeding_mode(bool mod) override {seeding_mode = mod;}
	void set_hits_map(std::map<unsigned int, SimpleHit3D>& hits_map) override{_hits_map = hits_map;}
	void set_verbose(int v) {verbose = v;}

	int run(std::vector<SimpleTrack3D>& output_tracks, std::vector<HelixKalmanState>& output_track_states, std::map<unsigned int, bool>& hits_used) override;	


 private:

	void set_detector_radii(std::vector<float>& radii) override;
	void set_detector_materials(std::vector<float>& materials);
	void set_input_tracks(std::vector<SimpleTrack3D>& input_tracks) override;
	void set_cylinder_kalman() override;

	int init() override;
	int process_tracks() override;
	int process_single_track(SimpleTrack3D& track) override;
	int process_single_triplet(SimpleTrack3D& track) override;
	int get_ca_tracks(std::vector<SimpleTrack3D>& output_tracks, std::vector<HelixKalmanState>& output_track_states) override;	

	int calculate_kappa_tangents(                        
			float x1, float y1, float z1, float x2, float y2, float z2, 
                        float x3, float y3, float z3,
                        float dx1, float dy1, float dz1, float dx2, float dy2, float dz2,
                        float dx3, float dy3, float dz3, 
                        float& kappa, float& dkappa,
                        float& ux_mid, float& uy_mid, float& ux_end, float& uy_end,
                        float& dzdl_1, float& dzdl_2, float& ddzdl_1, float& ddzdl_2) override;
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
                        float cur_chi2, float& chi2) override;	

	float kappa_to_pt(float kappa);
	float shift_phi_range(float phi) override;

	HelixHoughSpace* _hough_space;
	HelixKalmanFilter* _kalman;
	std::vector<SimpleTrack3D> in_tracks;
	std::vector<SimpleTrack3D> ca_tracks;
	std::vector<HelixKalmanState> ca_track_states;

	std::vector<unsigned int> temp_combo;
  	std::set<std::vector<unsigned int> > combos;
        std::vector<std::vector<SimpleHit3D> > layer_sorted;

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
	std::map<unsigned int, SimpleHit3D> _hits_map;

	double CAtime;
	double KALtime;	
	bool forward;
	bool remove_hits;
	bool remove_inner_hits;
	bool require_inner_hits;
	bool triplet_mode;
	bool seeding_mode;
	int verbose;

};

#endif
