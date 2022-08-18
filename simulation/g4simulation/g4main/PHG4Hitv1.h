// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4HITV1_H
#define G4MAIN_PHG4HITV1_H

#include "PHG4Hit.h"
#include "PHG4HitDefs.h"

#include <climits>  // for INT_MIN, ULONG_LONG_MAX
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>

class PHG4Hitv1 : public PHG4Hit
{
 public:
  PHG4Hitv1() = default;
  explicit PHG4Hitv1(const PHG4Hit* g4hit);
  ~PHG4Hitv1() override = default;
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;

  // The indices here represent the entry and exit points of the particle
  float get_x(const int i) const override { return x[i]; }
  float get_y(const int i) const override { return y[i]; }
  float get_z(const int i) const override { return z[i]; }
  float get_t(const int i) const override { return t[i]; }
  float get_edep() const override { return edep; }
  PHG4HitDefs::keytype get_hit_id() const override { return hitid; }
  int get_detid() const override;
  int get_shower_id() const override { return showerid; }
  int get_trkid() const override { return trackid; }

  void set_x(const int i, const float f) override { x[i] = f; }
  void set_y(const int i, const float f) override { y[i] = f; }
  void set_z(const int i, const float f) override { z[i] = f; }
  void set_t(const int i, const float f) override { t[i] = f; }
  void set_edep(const float f) override { edep = f; }
  void set_hit_id(const PHG4HitDefs::keytype i) override { hitid = i; }
  void set_shower_id(const int i) override { showerid = i; }
  void set_trkid(const int i) override { trackid = i; }

  void print() const override;

  bool has_property(const PROPERTY prop_id) const override;
  float get_property_float(const PROPERTY prop_id) const override;
  int get_property_int(const PROPERTY prop_id) const override;
  unsigned int get_property_uint(const PROPERTY prop_id) const override;
  void set_property(const PROPERTY prop_id, const float value) override;
  void set_property(const PROPERTY prop_id, const int value) override;
  void set_property(const PROPERTY prop_id, const unsigned int value) override;

  float get_px(const int i) const override;
  float get_py(const int i) const override;
  float get_pz(const int i) const override;
  float get_local_x(const int i) const override;
  float get_local_y(const int i) const override;
  float get_local_z(const int i) const override;
  float get_eion() const override { return get_property_float(prop_eion); }
  float get_light_yield() const override { return get_property_float(prop_light_yield); }
  float get_raw_light_yield() const override { return get_property_float(prop_raw_light_yield); }
  float get_path_length() const override { return get_property_float(prop_path_length); }
  unsigned int get_layer() const override { return get_property_uint(prop_layer); }
  int get_scint_id() const override { return get_property_int(prop_scint_id); }
  int get_row() const override { return get_property_int(prop_row); }
  int get_sector() const override { return get_property_int(prop_sector); }
  int get_strip_z_index() const override { return get_property_int(prop_strip_z_index); }
  int get_strip_y_index() const override { return get_property_int(prop_strip_y_index); }
  int get_ladder_z_index() const override { return get_property_int(prop_ladder_z_index); }
  int get_ladder_phi_index() const override { return get_property_int(prop_ladder_phi_index); }
  int get_index_i() const override { return get_property_int(prop_index_i); }
  int get_index_j() const override { return get_property_int(prop_index_j); }
  int get_index_k() const override { return get_property_int(prop_index_k); }
  int get_index_l() const override { return get_property_int(prop_index_l); }
  int get_hit_type() const override { return get_property_int(prop_hit_type); }

  void set_px(const int i, const float f) override;
  void set_py(const int i, const float f) override;
  void set_pz(const int i, const float f) override;
  void set_local_x(const int i, const float f) override;
  void set_local_y(const int i, const float f) override;
  void set_local_z(const int i, const float f) override;
  void set_eion(const float f) override { set_property(prop_eion, f); }
  void set_light_yield(const float f) override { set_property(prop_light_yield, f); }
  void set_raw_light_yield(const float f) override { set_property(prop_raw_light_yield, f); }
  void set_path_length(const float f) override { set_property(prop_path_length, f); }
  void set_layer(const unsigned int i) override { set_property(prop_layer, i); }
  void set_scint_id(const int i) override { set_property(prop_scint_id, i); }
  void set_row(const int i) override { set_property(prop_row, i); }
  void set_sector(const int i) override { set_property(prop_sector, i); }
  void set_strip_z_index(const int i) override { set_property(prop_strip_z_index, i); }
  void set_strip_y_index(const int i) override { set_property(prop_strip_y_index, i); }
  void set_ladder_z_index(const int i) override { set_property(prop_ladder_z_index, i); }
  void set_ladder_phi_index(const int i) override { set_property(prop_ladder_phi_index, i); }
  void set_index_i(const int i) override { set_property(prop_index_i, i); }
  void set_index_j(const int i) override { set_property(prop_index_j, i); }
  void set_index_k(const int i) override { set_property(prop_index_k, i); }
  void set_index_l(const int i) override { set_property(prop_index_l, i); }
  void set_hit_type(const int i) override { set_property(prop_hit_type, i); }

 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const override;
  void set_property_nocheck(const PROPERTY prop_id, const unsigned int ui) override { prop_map[prop_id] = ui; }
  // Store both the entry and exit points of the particle
  // Remember, particles do not always enter on the inner edge!
  float x[2] = {NAN, NAN};
  float y[2] = {NAN, NAN};
  float z[2] = {NAN, NAN};
  float t[2] = {NAN, NAN};
  PHG4HitDefs::keytype hitid = ULONG_LONG_MAX;
  int trackid = INT_MIN;
  int showerid = INT_MIN;
  float edep = NAN;

  //! storage types for additional property
  typedef uint8_t prop_id_t;
  typedef uint32_t prop_storage_t;
  typedef std::map<prop_id_t, prop_storage_t> prop_map_t;

  //! convert between 32bit inputs and storage type prop_storage_t
  union u_property
  {
    float fdata;
    int32_t idata;
    uint32_t uidata;

    u_property(int32_t in)
      : idata(in)
    {
    }
    u_property(uint32_t in)
      : uidata(in)
    {
    }
    u_property(float in)
      : fdata(in)
    {
    }
    u_property()
      : uidata(0)
    {
    }

    prop_storage_t get_storage() const { return uidata; }
  };

  //! container for additional property
  prop_map_t prop_map;

  ClassDefOverride(PHG4Hitv1, 2)
};

#endif
