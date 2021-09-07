#ifndef TRACKRECO_CELLULARAUTOMATON_H
#define TRACKRECO_CELLULARAUTOMATON_H


#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/HelixKalmanState.h>

#include <climits>
#include <cmath>
#include <map>
#include <set>

class HelixKalmanState;
class HelixHoughSpace;

class CellularAutomaton  {

 public:
  virtual ~CellularAutomaton() {};

  // The "standard PHObject response" functions...
  virtual void identify(std::ostream &os=std::cout) const {
    os << "CellularAutomaton base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int  isValid() const {return 0;}
  virtual CellularAutomaton* Clone() const {return nullptr;}

  virtual void set_hough_space(HelixHoughSpace*) {}
  virtual void set_mag_field(float) {}
  virtual void set_pt_rescale(float) {}
  virtual void set_n_layers(unsigned int) {}
  virtual void set_required_layers(unsigned int) {}
  virtual void set_ca_chi2(float) {}
  virtual void set_ca_chi2_layer(float) {}
  virtual void set_ca_phi_cut(float) {}
  virtual void set_ca_z_cut(float) {}
  virtual void set_ca_dcaxy_cut(float) {}
  virtual void set_propagate_forward(bool) {}
  virtual void set_remove_hits(bool) {}
  virtual void set_remove_inner_hits(bool) {}
  virtual void set_require_inner_hits(bool) {}
  virtual void set_triplet_mode(bool) {}
  virtual void set_seeding_mode(bool) {}
  virtual void set_hits_map(std::map<unsigned int, SimpleHit3D>&) {}

  virtual int run(std::vector<SimpleTrack3D>&, std::vector<HelixKalmanState>&, std::map<unsigned int, bool>&) {return 0;}

 private:

  virtual void set_detector_radii(std::vector<float>&) {}
  virtual void set_detector_material(std::vector<float>&) {}
  virtual void set_input_tracks(std::vector<SimpleTrack3D>&) {}
  virtual void set_cylinder_kalman() {}

  virtual int init() {return 0;}
  virtual int process_tracks() {return 0;}
  virtual int process_single_track(SimpleTrack3D&) {return 0;}
  virtual int process_single_triplet(SimpleTrack3D&) {return 0;}
  virtual int get_ca_tracks(std::vector<SimpleTrack3D>&, std::vector<HelixKalmanState>&) {return 0;}

  virtual int calculate_kappa_tangents(
                        float, float, float, float, float, float,
                        float, float, float,
                        float, float, float, float, float, float,
                        float, float, float,
                        float&, float&,
                        float&, float&, float&, float&,
                        float&, float&, float&, float&) {return 0;}

  virtual int calculate_kappa_tangents(
                        float, float, float, float, float, float,
                        float, float, float,
                        float, float, float, float, float, float,
                        float, float, float,
                        float&, float&,
                        float&, float&, float&, float&,
                        float&, float&, float&, float&,
                        float, float,
                        float, float, float, float,
                        float, float&) {return 0;}

  virtual float shift_phi_range(float){return 0;};


 protected:
  CellularAutomaton(){}

};

#endif
