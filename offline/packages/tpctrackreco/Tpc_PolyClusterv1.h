#pragma once

#include "Tpc_PolyCluster.h"

#include <iostream>
#include <vector>

class Tpc_PolyClusterv1 : public Tpc_PolyCluster
{
 public:
  Tpc_PolyClusterv1();
  ~Tpc_PolyClusterv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new Tpc_PolyClusterv1(*this); }

  unsigned int get_event() const override { return m_event; }
  unsigned int get_cluster_id() const override { return m_cluster_id; }
  unsigned int get_source_assembled_track_id() const override { return m_source_assembled_track_id; }
  int get_side() const override { return m_side; }
  unsigned int get_nhits() const override { return static_cast<unsigned int>(m_hit_indices.size()); }

  double get_centroid_x() const override { return m_centroid_x; }
  double get_centroid_y() const override { return m_centroid_y; }
  double get_centroid_z() const override { return m_centroid_z; }
  double get_rms_x() const override { return m_rms_x; }
  double get_rms_y() const override { return m_rms_y; }
  double get_rms_z() const override { return m_rms_z; }
  double get_adc() const override { return m_adc; }
  unsigned int get_phi_width() const override { return m_phi_width; }
  unsigned int get_time_width() const override { return m_time_width; }
  double get_phase() const override { return m_phase; }

  void set_event(unsigned int v) override { m_event = v; }
  void set_cluster_id(unsigned int v) override { m_cluster_id = v; }
  void set_source_assembled_track_id(unsigned int v) override { m_source_assembled_track_id = v; }
  void set_side(int v) override { m_side = v; }
  void set_centroid_x(double v) override { m_centroid_x = v; }
  void set_centroid_y(double v) override { m_centroid_y = v; }
  void set_centroid_z(double v) override { m_centroid_z = v; }
  void set_rms_x(double v) override { m_rms_x = v; }
  void set_rms_y(double v) override { m_rms_y = v; }
  void set_rms_z(double v) override { m_rms_z = v; }
  void set_adc(double v) override { m_adc = v; }
  void set_phi_width(unsigned int v) override { m_phi_width = v; }
  void set_time_width(unsigned int v) override { m_time_width = v; }
  void set_phase(double v) override { m_phase = v; }

  void add_hit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk,
               double x, double y, double z) override
  {
    m_hit_indices.emplace_back(hsk, hk);
    m_hit_x.push_back(x);
    m_hit_y.push_back(y);
    m_hit_z.push_back(z);
  }

  unsigned int size_hits() const override { return static_cast<unsigned int>(m_hit_indices.size()); }
  HitIndex get_hit_index(unsigned int i) const override
  {
    if (i >= m_hit_indices.size()) return {0, 0};
    return m_hit_indices[i];
  }
  double get_hit_x(unsigned int i) const override { return i < m_hit_x.size() ? m_hit_x[i] : 0.0; }
  double get_hit_y(unsigned int i) const override { return i < m_hit_y.size() ? m_hit_y[i] : 0.0; }
  double get_hit_z(unsigned int i) const override { return i < m_hit_z.size() ? m_hit_z[i] : 0.0; }
  const std::vector<HitIndex>& get_hit_indices() const override { return m_hit_indices; }

 private:
  unsigned int m_event {0};
  unsigned int m_cluster_id {0};
  unsigned int m_source_assembled_track_id {0};
  int m_side {0};

  double m_centroid_x {0.0};
  double m_centroid_y {0.0};
  double m_centroid_z {0.0};
  double m_rms_x {0.0};
  double m_rms_y {0.0};
  double m_rms_z {0.0};
  double m_adc {0.0};
  unsigned int m_phi_width {0};
  unsigned int m_time_width {0};
  double m_phase {0.0};

  std::vector<HitIndex> m_hit_indices;
  std::vector<double> m_hit_x;
  std::vector<double> m_hit_y;
  std::vector<double> m_hit_z;

  ClassDefOverride(Tpc_PolyClusterv1, 1)
};
