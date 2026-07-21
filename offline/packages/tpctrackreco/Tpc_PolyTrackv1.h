#pragma once

#include "Tpc_PolyTrack.h"

#include <iostream>
#include <vector>

class Tpc_PolyTrackv1 : public Tpc_PolyTrack
{
 public:
  Tpc_PolyTrackv1();
  ~Tpc_PolyTrackv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new Tpc_PolyTrackv1(*this); }

  unsigned int get_event() const override { return m_event; }
  unsigned int get_track_id() const override { return m_track_id; }
  unsigned int get_source_assembled_track_id() const override { return m_source_assembled_track_id; }
  int get_fit_status() const override { return m_fit_status; }
  unsigned int get_nclusters() const override { return m_nclusters; }
  double get_x() const override { return m_x; }
  double get_y() const override { return m_y; }
  double get_z() const override { return m_z; }
  double get_px() const override { return m_px; }
  double get_py() const override { return m_py; }
  double get_pz() const override { return m_pz; }
  double get_charge() const override { return m_charge; }
  double get_chi2() const override { return m_chi2; }
  double get_ndf() const override { return m_ndf; }
  double get_dedx() const override { return m_dedx; }
  double get_cov(unsigned int i, unsigned int j) const override
  {
    const unsigned int idx = 6 * i + j;
    if (i >= 6 || j >= 6 || idx >= m_cov.size()) return 0.0;
    return m_cov[idx];
  }

  void set_event(unsigned int v) override { m_event = v; }
  void set_track_id(unsigned int v) override { m_track_id = v; }
  void set_source_assembled_track_id(unsigned int v) override { m_source_assembled_track_id = v; }
  void set_fit_status(int v) override { m_fit_status = v; }
  void set_nclusters(unsigned int v) override { m_nclusters = v; }
  void set_x(double v) override { m_x = v; }
  void set_y(double v) override { m_y = v; }
  void set_z(double v) override { m_z = v; }
  void set_px(double v) override { m_px = v; }
  void set_py(double v) override { m_py = v; }
  void set_pz(double v) override { m_pz = v; }
  void set_charge(double v) override { m_charge = v; }
  void set_chi2(double v) override { m_chi2 = v; }
  void set_ndf(double v) override { m_ndf = v; }
  void set_dedx(double v) override { m_dedx = v; }
  void set_cov(unsigned int i, unsigned int j, double v) override
  {
    if (i >= 6 || j >= 6) return;
    if (m_cov.size() != 36) m_cov.assign(36, 0.0);
    m_cov[6 * i + j] = v;
    m_cov[6 * j + i] = v;
  }

 private:
  unsigned int m_event {0};
  unsigned int m_track_id {0};
  unsigned int m_source_assembled_track_id {0};
  int m_fit_status {0};
  unsigned int m_nclusters {0};
  double m_x {0.0};
  double m_y {0.0};
  double m_z {0.0};
  double m_px {0.0};
  double m_py {0.0};
  double m_pz {0.0};
  double m_charge {0.0};
  double m_chi2 {0.0};
  double m_ndf {0.0};
  double m_dedx {0.0};
  std::vector<double> m_cov;

  ClassDefOverride(Tpc_PolyTrackv1, 2)
};
