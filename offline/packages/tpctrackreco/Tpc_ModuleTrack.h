#pragma once

#include <trackbase/TrkrDefs.h>

#include <phool/PHObject.h>

#include <iostream>
#include <utility>
#include <vector>

// A single in-module track.
// It stores only identity/topology and (hitsetkey, hitkey) references into
// TrkrHitSetContainer.  No fit result is stored here; run a separate Tpc_FittingTools
// module if fit parameters are needed downstream.
class Tpc_ModuleTrack : public PHObject
{
 public:
  Tpc_ModuleTrack() = default;
  ~Tpc_ModuleTrack() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_ModuleTrack base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  // --- Identity ---
  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_track_id() const { return 0; }

  virtual unsigned int get_region() const { return 0; }
  virtual unsigned int get_sector() const { return 0; }
  virtual int get_side() const { return 0; }

  // --- Topology ---
  virtual unsigned int get_nblobs() const { return 0; }
  virtual unsigned int get_nrawhits() const { return 0; }

  virtual unsigned int get_first_layer() const { return 0; }
  virtual unsigned int get_last_layer() const { return 0; }

  // --- Setters ---
  virtual void set_event(unsigned int) {}
  virtual void set_track_id(unsigned int) {}

  virtual void set_region(unsigned int) {}
  virtual void set_sector(unsigned int) {}
  virtual void set_side(int) {}

  virtual void set_nblobs(unsigned int) {}
  virtual void set_nrawhits(unsigned int) {}

  virtual void set_first_layer(unsigned int) {}
  virtual void set_last_layer(unsigned int) {}

  // --- Hit index access ---
  using HitIndex = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;

  virtual void add_hit_index(TrkrDefs::hitsetkey /*hsk*/, TrkrDefs::hitkey /*hk*/) {}

  virtual unsigned int size_hit_indices() const { return 0; }

  virtual HitIndex get_hit_index(unsigned int /*i*/) const
  {
    return {0, 0};
  }

  virtual const std::vector<HitIndex>& get_hit_indices() const
  {
    static const std::vector<HitIndex> empty_indices;
    return empty_indices;
  }

 private:
  ClassDefOverride(Tpc_ModuleTrack, 1)
};
