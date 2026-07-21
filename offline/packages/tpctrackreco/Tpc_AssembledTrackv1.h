#pragma once

#include "Tpc_AssembledTrack.h"

#include <iostream>
#include <vector>

class Tpc_AssembledTrackv1 : public Tpc_AssembledTrack
{
 public:
  Tpc_AssembledTrackv1();
  ~Tpc_AssembledTrackv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new Tpc_AssembledTrackv1(*this); }

  unsigned int get_event() const override { return m_event; }
  unsigned int get_track_id() const override { return m_track_id; }
  int get_side() const override { return m_side; }

  unsigned int get_nsegments() const override { return m_nsegments; }
  unsigned int get_nblobs() const override { return m_nblobs; }
  unsigned int get_nrawhits() const override { return m_nrawhits; }
  unsigned int get_first_layer() const override { return m_first_layer; }
  unsigned int get_last_layer() const override { return m_last_layer; }
  unsigned int get_first_sector() const override { return m_first_sector; }
  unsigned int get_last_sector() const override { return m_last_sector; }
  unsigned int get_first_region() const override { return m_first_region; }
  unsigned int get_last_region() const override { return m_last_region; }

  double get_phi_slope() const override { return m_phi_slope; }
  double get_phi_intercept() const override { return m_phi_intercept; }
  double get_tbin_slope() const override { return m_tbin_slope; }
  double get_tbin_intercept() const override { return m_tbin_intercept; }
  double get_chi2_phi() const override { return m_chi2_phi; }
  double get_chi2_tbin() const override { return m_chi2_tbin; }
  int get_ndof_phi() const override { return m_ndof_phi; }
  int get_ndof_tbin() const override { return m_ndof_tbin; }

  int get_vertex_valid() const override { return m_vertex_valid; }
  double get_vertex_x() const override { return m_vertex_x; }
  double get_vertex_y() const override { return m_vertex_y; }
  double get_vertex_r() const override { return m_vertex_r; }
  double get_vertex_phi() const override { return m_vertex_phi; }
  double get_vertex_tbin() const override { return m_vertex_tbin; }
  unsigned int get_vertex_npairs() const override { return m_vertex_npairs; }
  double get_vertex_quality() const override { return m_vertex_quality; }

  int get_seed_valid() const override { return m_seed_valid; }
  double get_seed_x() const override { return m_seed_x; }
  double get_seed_y() const override { return m_seed_y; }
  double get_seed_z() const override { return m_seed_z; }
  double get_seed_px() const override { return m_seed_px; }
  double get_seed_py() const override { return m_seed_py; }
  double get_seed_pz() const override { return m_seed_pz; }
  double get_seed_cov(unsigned int i, unsigned int j) const override
  {
    const unsigned int idx = 6 * i + j;
    if (i >= 6 || j >= 6 || idx >= m_seed_cov.size()) return 0.0;
    return m_seed_cov[idx];
  }

  void set_vertex_valid(int v) override { m_vertex_valid = v; }
  void set_vertex_x(double v) override { m_vertex_x = v; }
  void set_vertex_y(double v) override { m_vertex_y = v; }
  void set_vertex_r(double v) override { m_vertex_r = v; }
  void set_vertex_phi(double v) override { m_vertex_phi = v; }
  void set_vertex_tbin(double v) override { m_vertex_tbin = v; }
  void set_vertex_npairs(unsigned int v) override { m_vertex_npairs = v; }
  void set_vertex_quality(double v) override { m_vertex_quality = v; }

  void set_seed_valid(int v) override { m_seed_valid = v; }
  void set_seed_x(double v) override { m_seed_x = v; }
  void set_seed_y(double v) override { m_seed_y = v; }
  void set_seed_z(double v) override { m_seed_z = v; }
  void set_seed_px(double v) override { m_seed_px = v; }
  void set_seed_py(double v) override { m_seed_py = v; }
  void set_seed_pz(double v) override { m_seed_pz = v; }
  void set_seed_cov(unsigned int i, unsigned int j, double v) override
  {
    if (i >= 6 || j >= 6) return;
    if (m_seed_cov.size() != 36) m_seed_cov.assign(36, 0.0);
    m_seed_cov[6 * i + j] = v;
    m_seed_cov[6 * j + i] = v;
  }

  void set_event(unsigned int v) override { m_event = v; }
  void set_track_id(unsigned int v) override { m_track_id = v; }
  void set_side(int v) override { m_side = v; }
  void set_nsegments(unsigned int v) override { m_nsegments = v; }
  void set_nblobs(unsigned int v) override { m_nblobs = v; }
  void set_nrawhits(unsigned int v) override { m_nrawhits = v; }
  void set_first_layer(unsigned int v) override { m_first_layer = v; }
  void set_last_layer(unsigned int v) override { m_last_layer = v; }
  void set_first_sector(unsigned int v) override { m_first_sector = v; }
  void set_last_sector(unsigned int v) override { m_last_sector = v; }
  void set_first_region(unsigned int v) override { m_first_region = v; }
  void set_last_region(unsigned int v) override { m_last_region = v; }
  void set_phi_slope(double v) override { m_phi_slope = v; }
  void set_phi_intercept(double v) override { m_phi_intercept = v; }
  void set_tbin_slope(double v) override { m_tbin_slope = v; }
  void set_tbin_intercept(double v) override { m_tbin_intercept = v; }
  void set_chi2_phi(double v) override { m_chi2_phi = v; }
  void set_chi2_tbin(double v) override { m_chi2_tbin = v; }
  void set_ndof_phi(int v) override { m_ndof_phi = v; }
  void set_ndof_tbin(int v) override { m_ndof_tbin = v; }

  void add_source_track(unsigned int track_id,
                        unsigned int region,
                        unsigned int sector) override
  {
    m_source_track_ids.push_back(track_id);
    m_source_regions.push_back(region);
    m_source_sectors.push_back(sector);
  }

  unsigned int size_source_tracks() const override
  {
    return static_cast<unsigned int>(m_source_track_ids.size());
  }

  unsigned int get_source_track_id(unsigned int i) const override
  {
    if (i >= m_source_track_ids.size()) return 0;
    return m_source_track_ids[i];
  }

  unsigned int get_source_region(unsigned int i) const override
  {
    if (i >= m_source_regions.size()) return 0;
    return m_source_regions[i];
  }

  unsigned int get_source_sector(unsigned int i) const override
  {
    if (i >= m_source_sectors.size()) return 0;
    return m_source_sectors[i];
  }

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
  int m_side {0};

  unsigned int m_nsegments {0};
  unsigned int m_nblobs {0};
  unsigned int m_nrawhits {0};
  unsigned int m_first_layer {0};
  unsigned int m_last_layer {0};
  unsigned int m_first_sector {0};
  unsigned int m_last_sector {0};
  unsigned int m_first_region {0};
  unsigned int m_last_region {0};

  double m_phi_slope {0.0};
  double m_phi_intercept {0.0};
  double m_tbin_slope {0.0};
  double m_tbin_intercept {0.0};
  double m_chi2_phi {0.0};
  double m_chi2_tbin {0.0};
  int m_ndof_phi {0};
  int m_ndof_tbin {0};

  int m_vertex_valid {0};
  double m_vertex_x {0.0};
  double m_vertex_y {0.0};
  double m_vertex_r {0.0};
  double m_vertex_phi {0.0};
  double m_vertex_tbin {0.0};
  unsigned int m_vertex_npairs {0};
  double m_vertex_quality {0.0};

  int m_seed_valid {0};
  double m_seed_x {0.0};
  double m_seed_y {0.0};
  double m_seed_z {0.0};
  double m_seed_px {0.0};
  double m_seed_py {0.0};
  double m_seed_pz {0.0};
  std::vector<double> m_seed_cov;

  std::vector<unsigned int> m_source_track_ids;
  std::vector<unsigned int> m_source_regions;
  std::vector<unsigned int> m_source_sectors;
  std::vector<HitIndex> m_hit_indices;

  ClassDefOverride(Tpc_AssembledTrackv1, 2)
};
