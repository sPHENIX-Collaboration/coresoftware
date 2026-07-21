#pragma once

#include "Tpc_ModuleTrack.h"

#include <iostream>
#include <utility>
#include <vector>

class Tpc_ModuleTrackv1 : public Tpc_ModuleTrack
{
 public:
  Tpc_ModuleTrackv1();
  ~Tpc_ModuleTrackv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new Tpc_ModuleTrackv1(*this); }

  // --- Identity ---
  unsigned int get_event() const override { return m_event; }
  unsigned int get_track_id() const override { return m_track_id; }

  unsigned int get_region() const override { return m_region; }
  unsigned int get_sector() const override { return m_sector; }
  int get_side() const override { return m_side; }

  // --- Topology ---
  unsigned int get_nblobs() const override { return m_nblobs; }
  unsigned int get_nrawhits() const override { return m_nrawhits; }

  unsigned int get_first_layer() const override { return m_first_layer; }
  unsigned int get_last_layer() const override { return m_last_layer; }

  // --- Setters ---
  void set_event(unsigned int v) override { m_event = v; }
  void set_track_id(unsigned int v) override { m_track_id = v; }

  void set_region(unsigned int v) override { m_region = v; }
  void set_sector(unsigned int v) override { m_sector = v; }
  void set_side(int v) override { m_side = v; }

  void set_nblobs(unsigned int v) override { m_nblobs = v; }
  void set_nrawhits(unsigned int v) override { m_nrawhits = v; }

  void set_first_layer(unsigned int v) override { m_first_layer = v; }
  void set_last_layer(unsigned int v) override { m_last_layer = v; }

  // --- Hit index access ---
  void add_hit_index(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk) override
  {
    m_hit_indices.emplace_back(hsk, hk);
  }

  unsigned int size_hit_indices() const override
  {
    return static_cast<unsigned int>(m_hit_indices.size());
  }

  HitIndex get_hit_index(unsigned int i) const override
  {
    if (i >= m_hit_indices.size()) return {0, 0};
    return m_hit_indices[i];
  }

  const std::vector<HitIndex>& get_hit_indices() const override
  {
    return m_hit_indices;
  }

 private:
  unsigned int m_event {0};
  unsigned int m_track_id {0};

  unsigned int m_region {0};
  unsigned int m_sector {0};
  int m_side {0};

  unsigned int m_nblobs {0};
  unsigned int m_nrawhits {0};

  unsigned int m_first_layer {0};
  unsigned int m_last_layer {0};

  std::vector<HitIndex> m_hit_indices;

  ClassDefOverride(Tpc_ModuleTrackv1, 3)
};
