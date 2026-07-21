#pragma once

#include "Tpc_PolyTrackVertex.h"

#include <iostream>

class Tpc_PolyTrackVertexv1 : public Tpc_PolyTrackVertex
{
 public:
  Tpc_PolyTrackVertexv1();
  ~Tpc_PolyTrackVertexv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new Tpc_PolyTrackVertexv1(*this); }

  unsigned int get_track_id() const override { return m_track_id; }
  unsigned int get_source_assembled_track_id() const override { return m_source_assembled_track_id; }
  double get_dca2d() const override { return m_dca2d; }
  double get_z0() const override { return m_z0; }
  int get_pca_valid() const override { return m_pca_valid; }
  double get_pca_x() const override { return m_pca_x; }
  double get_pca_y() const override { return m_pca_y; }
  double get_pca_z() const override { return m_pca_z; }
  double get_pca_radius() const override { return m_pca_radius; }
  double get_pca_phi() const override { return m_pca_phi; }

  void set_track_id(unsigned int v) override { m_track_id = v; }
  void set_source_assembled_track_id(unsigned int v) override { m_source_assembled_track_id = v; }
  void set_dca2d(double v) override { m_dca2d = v; }
  void set_z0(double v) override { m_z0 = v; }
  void set_pca_valid(int v) override { m_pca_valid = v; }
  void set_pca_x(double v) override { m_pca_x = v; }
  void set_pca_y(double v) override { m_pca_y = v; }
  void set_pca_z(double v) override { m_pca_z = v; }
  void set_pca_radius(double v) override { m_pca_radius = v; }
  void set_pca_phi(double v) override { m_pca_phi = v; }

 private:
  unsigned int m_track_id {0};
  unsigned int m_source_assembled_track_id {0};
  double m_dca2d {0.0};
  double m_z0 {0.0};
  int m_pca_valid {0};
  double m_pca_x {0.0};
  double m_pca_y {0.0};
  double m_pca_z {0.0};
  double m_pca_radius {0.0};
  double m_pca_phi {0.0};

  ClassDefOverride(Tpc_PolyTrackVertexv1, 1)
};
