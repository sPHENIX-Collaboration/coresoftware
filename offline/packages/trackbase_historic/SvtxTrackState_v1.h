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
  ~SvtxTrackState_v1() override {}

  // The "standard PHObject response" functions...
  void identify(std::ostream &os = std::cout) const override;
  void Reset() override { *this = SvtxTrackState_v1(0.0); }
  int isValid() const override { return 1; }
  PHObject *CloneMe() const override { return new SvtxTrackState_v1(*this); }

  float get_pathlength() const override { return _pathlength; }

  float get_x() const override { return _pos[0]; }
  void set_x(float x) override { _pos[0] = x; }

  float get_y() const override { return _pos[1]; }
  void set_y(float y) override { _pos[1] = y; }

  float get_z() const override { return _pos[2]; }
  void set_z(float z) override { _pos[2] = z; }

  float get_pos(unsigned int i) const override { return _pos[i]; }

  float get_px() const override { return _mom[0]; }
  void set_px(float px) override { _mom[0] = px; }

  float get_py() const override { return _mom[1]; }
  void set_py(float py) override { _mom[1] = py; }

  float get_pz() const override { return _mom[2]; }
  void set_pz(float pz) override { _mom[2] = pz; }

  float get_mom(unsigned int i) const override { return _mom[i]; }

  float get_p() const override { return sqrt(pow(get_px(), 2) + pow(get_py(), 2) + pow(get_pz(), 2)); }
  float get_pt() const override { return sqrt(pow(get_px(), 2) + pow(get_py(), 2)); }
  float get_eta() const override { return asinh(get_pz() / get_pt()); }
  float get_phi() const override { return atan2(get_py(), get_px()); }

  float get_error(unsigned int i, unsigned int j) const override;
  void set_error(unsigned int i, unsigned int j, float value) override;

  std::string get_name() const override { return state_name; }
  void set_name(const std::string &name) override { state_name = name; }

  float get_rphi_error() const override;
  float get_phi_error() const override;
  float get_z_error() const override;

  //@}

 private:
  float _pathlength;
  float _pos[3];
  float _mom[3];
  float _covar[21];  //  6x6 triangular packed storage

  std::string state_name;

  ClassDefOverride(SvtxTrackState_v1, 1)
};

#endif
