// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CELLV1_H
#define G4DETECTORS_PHG4CELLV1_H

#include "PHG4Cell.h"
#include "PHG4CellDefs.h"

#include <g4main/PHG4HitDefs.h>  // for keytype

#include <cstdint>
#include <iostream>
#include <map>
#include <type_traits>  // for __decay_and_strip<>::__type
#include <utility>      // for make_pair

class PHG4Cellv1 : public PHG4Cell
{
 public:
  PHG4Cellv1() {}
  explicit PHG4Cellv1(const PHG4CellDefs::keytype g4cellid);
  ~PHG4Cellv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;

  void set_cellid(const PHG4CellDefs::keytype i) override { cellid = i; }

  PHG4CellDefs::keytype get_cellid() const override { return cellid; }
  bool has_binning(const PHG4CellDefs::CellBinning binning) const override;
  short int get_detid() const override;

  using PHG4Cell::add_edep;

  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep) override;
  void add_shower_edep(const int g4showerid, const float edep) override;

  EdepConstRange get_g4hits() override
  {
    return std::make_pair(hitedeps.begin(), hitedeps.end());
  }

  ShowerEdepConstRange get_g4showers() override
  {
    return std::make_pair(showeredeps.begin(), showeredeps.end());
  }

  void add_edep(const float f) override { add_property(prop_edep, f); }
  double get_edep() const override { return get_property_float(prop_edep); }

  void add_eion(const float f) override { add_property(prop_eion, f); }
  double get_eion() const override { return get_property_float(prop_eion); }

  void add_light_yield(const float f) override { add_property(prop_light_yield, f); }
  float get_light_yield() const override { return get_property_float(prop_light_yield); }

  void add_raw_light_yield(const float f) override { add_property(prop_raw_light_yield, f); }
  float get_raw_light_yield() const override { return get_property_float(prop_raw_light_yield); }

  void set_chip_index(const int i) override { set_property(prop_chip_index, i); }
  int get_chip_index() const override { return get_property_int(prop_chip_index); }

  void set_half_stave_index(const int i) override { set_property(prop_half_stave_index, i); }
  int get_half_stave_index() const override { return get_property_int(prop_half_stave_index); }

  void set_ladder_phi_index(const int i) override { set_property(prop_ladder_phi_index, i); }
  int get_ladder_phi_index() const override { return get_property_int(prop_ladder_phi_index); }

  void set_ladder_z_index(const int i) override { set_property(prop_ladder_z_index, i); }
  int get_ladder_z_index() const override { return get_property_int(prop_ladder_z_index); }

  void set_module_index(const int i) override { set_property(prop_module_index, i); }
  int get_module_index() const override { return get_property_int(prop_module_index); }

  void set_phibin(const int i) override { set_property(prop_phibin, i); }
  int get_phibin() const override { return get_property_int(prop_phibin); }

  void set_pixel_index(const int i) override { set_property(prop_pixel_index, i); }
  int get_pixel_index() const override { return get_property_int(prop_pixel_index); }

  void set_stave_index(const int i) override { set_property(prop_stave_index, i); }
  int get_stave_index() const override { return get_property_int(prop_stave_index); }

  void set_zbin(const int i) override { set_property(prop_zbin, i); }
  int get_zbin() const override { return get_property_int(prop_zbin); }

  void print() const override;

  bool has_property(const PROPERTY prop_id) const override;
  float get_property_float(const PROPERTY prop_id) const override;
  int get_property_int(const PROPERTY prop_id) const override;
  unsigned int get_property_uint(const PROPERTY prop_id) const override;
  void add_property(const PROPERTY prop_id, const float value);
  void add_property(const PROPERTY prop_id, const int value);
  void add_property(const PROPERTY prop_id, const unsigned int value);
  void set_property(const PROPERTY prop_id, const float value) override;
  void set_property(const PROPERTY prop_id, const int value) override;
  void set_property(const PROPERTY prop_id, const unsigned int value) override;

 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const override;
  void set_property_nocheck(const PROPERTY prop_id, const unsigned int ui) override { prop_map[prop_id] = ui; }

  PHG4CellDefs::keytype cellid = ~0x0;
  EdepMap hitedeps;
  ShowerEdepMap showeredeps;

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

  ClassDefOverride(PHG4Cellv1, 3)
};

#endif
