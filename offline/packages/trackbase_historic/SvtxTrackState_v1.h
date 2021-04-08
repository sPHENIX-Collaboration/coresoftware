#ifndef TRACKBASEHISTORIC_SVTXTRACKSTATEV1_H
#define TRACKBASEHISTORIC_SVTXTRACKSTATEV1_H

#include "SvtxTrackState.h"

#include <cmath>
#include <iostream>
#include <string>  // for string, basic_string

class PHObject;

class SvtxTrackState_v1 : public SvtxTrackState
{
 public:
  SvtxTrackState_v1(float pathlength = 0.0);
  virtual ~SvtxTrackState_v1() {}

  // The "standard PHObject response" functions...
  void identify(std::ostream &os = std::cout) const;
  void Reset() { *this = SvtxTrackState_v1(0.0); }
  int isValid() const { return 1; }
  PHObject *CloneMe() const { return new SvtxTrackState_v1(*this); }

  float get_pathlength() const { return _pathlength; }

  float get_x() const { return _pos[0]; }
  void set_x(float x) { _pos[0] = x; }

  float get_y() const { return _pos[1]; }
  void set_y(float y) { _pos[1] = y; }

  float get_z() const { return _pos[2]; }
  void set_z(float z) { _pos[2] = z; }

  float get_pos(unsigned int i) const { return _pos[i]; }

  float get_px() const { return _mom[0]; }
  void set_px(float px) { _mom[0] = px; }

  float get_py() const { return _mom[1]; }
  void set_py(float py) { _mom[1] = py; }

  float get_pz() const { return _mom[2]; }
  void set_pz(float pz) { _mom[2] = pz; }

  float get_mom(unsigned int i) const { return _mom[i]; }

  float get_p() const { return sqrt(pow(get_px(), 2) + pow(get_py(), 2) + pow(get_pz(), 2)); }
  float get_pt() const { return sqrt(pow(get_px(), 2) + pow(get_py(), 2)); }
  float get_eta() const { return asinh(get_pz() / get_pt()); }
  float get_phi() const { return atan2(get_py(), get_px()); }

  float get_error(unsigned int i, unsigned int j) const;
  void set_error(unsigned int i, unsigned int j, float value);

  std::string get_name() const { return state_name; }
  void set_name(const std::string &name) { state_name = name; }

  virtual float get_rphi_error() const;
  virtual float get_phi_error() const;
  virtual float get_z_error() const;

  //@}

 private:
  float _pathlength;
  float _pos[3];
  float _mom[3];
  float _covar[21];  //  6x6 triangular packed storage

  std::string state_name;

  ClassDef(SvtxTrackState_v1, 1)
};

#endif
