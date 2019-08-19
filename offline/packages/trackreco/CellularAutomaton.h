#ifndef TRACKRECO_CELLULARAUTOMATON_H
#define TRACKRECO_CELLULARAUTOMATON_H

#if !defined(__CINT__) || defined(__CLING__)
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/HelixKalmanState.h>
#endif

#include "Track3D.h"
#include "HelixHoughSpace.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <climits>
#include <cmath>
#include <map>
#include <set>

class Track3D;
class HelixKalmanState;
class HelixHoughSpace;

class CellularAutomaton : public PHObject {

 public:
  virtual ~CellularAutomaton() {};

  // The "standard PHObject response" functions...
  virtual void identify(std::ostream &os=std::cout) const {
    os << "CellularAutomaton base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int  isValid() const {return 0;}
  virtual CellularAutomaton* Clone() const {return nullptr;}

  virtual void set_hough_space(HelixHoughSpace* hough_space) {}
  virtual void set_mag_field(float mag_field) {}
  virtual void set_pt_rescale(float pt_rescale) {}
  virtual void set_n_layers(unsigned int n_layers) {}
  virtual void set_required_layers(unsigned int req_layers) {}
  virtual void set_ca_chi2(float ca_chi2) {}
  virtual void set_ca_chi2_layer(float chi2_layer_cut) {}
  virtual void set_ca_phi_cut(float phi_cut) {}
  virtual void set_ca_z_cut(float z_cut) {}
  virtual void set_ca_dcaxy_cut(float dcaxy_cut) {}
  virtual void set_propagate_forward(bool fwd) {}
  virtual void set_remove_hits(bool remove) {}
  virtual void set_remove_inner_hits(bool remove_inner) {}
  virtual void set_require_inner_hits(bool require_inner) {}
  virtual void set_triplet_mode(bool mod) {}
  virtual void set_seeding_mode(bool mod) {}
  virtual void set_hits_map(std::map<unsigned int, SimpleHit3D>& hits_map) {}

  virtual int run(std::vector<Track3D>& output_tracks, std::vector<HelixKalmanState>& output_track_states, std::map<unsigned int, bool>& hits_used) {return 0;}

 private:

  virtual void set_detector_radii(std::vector<float>& radius) {}
  virtual void set_detector_material(std::vector<float>& material) {}
  virtual void set_input_tracks(std::vector<Track3D>& input_tracks) {}
  virtual void set_cylinder_kalman() {}

  virtual int init() {return 0;}
  virtual int process_tracks() {return 0;}
  virtual int process_single_track(Track3D& track) {return 0;}
  virtual int process_single_triplet(Track3D& track) {return 0;}
  virtual int get_ca_tracks(std::vector<Track3D>& output_tracks, std::vector<HelixKalmanState>& output_track_states) {return 0;}

  virtual int calculate_kappa_tangents(
                        float x1, float y1, float z1, float x2, float y2, float z2,
                        float x3, float y3, float z3,
                        float dx1, float dy1, float dz1, float dx2, float dy2, float dz2,
                        float dx3, float dy3, float dz3,
                        float& kappa, float& dkappa,
                        float& ux_mid, float& uy_mid, float& ux_end, float&uy_end,
                        float& dzdl_1, float& dzdl_2, float& ddzdl_1, float& ddzdl_2) {return 0;}

  virtual int calculate_kappa_tangents(
                        float x1, float y1, float z1, float x2, float y2, float z2,
                        float x3, float y3, float z3,
                        float dx1, float dy1, float dz1, float dx2, float dy2, float dz2,
                        float dx3, float dy3, float dz3,
                        float& kappa, float& dkappa,
                        float& ux_mid, float& uy_mid, float& ux_end, float& uy_end,
                        float& dzdl_1, float& dzdl_2, float& ddzdl_1, float& ddzdl_2,
                        float ca_sin_ang_cut, float inv_cos_ang_diff,
                        float cur_kappa, float cur_dkappa, float cur_ux, float cur_uy,
                        float cur_chi2, float& chi2) {return 0;}

  virtual float shift_phi_range(float phi){return 0;};


 protected:
  CellularAutomaton(){}

  ClassDef(CellularAutomaton,1);
};

#endif
