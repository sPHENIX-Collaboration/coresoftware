#ifndef TPCTRACKRECO_TPCCROSSINGDECISIONV1_H
#define TPCTRACKRECO_TPCCROSSINGDECISIONV1_H


#include "TpcCrossingDecision.h"

#include <iostream>
#include <limits>

class TpcCrossingDecisionv1 : public TpcCrossingDecision
{
 public:
  TpcCrossingDecisionv1();
  ~TpcCrossingDecisionv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new TpcCrossingDecisionv1(*this); }

  unsigned int get_assembled_track_id() const override { return m_assembled_track_id; }
  void set_assembled_track_id(unsigned int v) override { m_assembled_track_id = v; }

  short get_selected_crossing() const override { return m_selected_crossing; }
  void set_selected_crossing(short v) override { m_selected_crossing = v; }

  unsigned int get_silicon_vertex_id() const override { return m_silicon_vertex_id; }
  void set_silicon_vertex_id(unsigned int v) override { m_silicon_vertex_id = v; }

  float get_tpc_z0() const override { return m_tpc_z0; }
  void set_tpc_z0(float v) override { m_tpc_z0 = v; }

  float get_silicon_vertex_z() const override { return m_silicon_vertex_z; }
  void set_silicon_vertex_z(float v) override { m_silicon_vertex_z = v; }

  float get_delta_z() const override { return m_delta_z; }
  void set_delta_z(float v) override { m_delta_z = v; }

  float get_best_abs_delta_z() const override { return m_best_abs_delta_z; }
  void set_best_abs_delta_z(float v) override { m_best_abs_delta_z = v; }

  float get_second_best_abs_delta_z() const override { return m_second_best_abs_delta_z; }
  void set_second_best_abs_delta_z(float v) override { m_second_best_abs_delta_z = v; }

  unsigned char get_selected_tier() const override { return m_selected_tier; }
  void set_selected_tier(unsigned char v) override { m_selected_tier = v; }

  float get_selected_score() const override { return m_selected_score; }
  void set_selected_score(float v) override { m_selected_score = v; }

  unsigned short get_number_of_available_crossings() const override { return m_number_of_available_crossings; }
  void set_number_of_available_crossings(unsigned short v) override { m_number_of_available_crossings = v; }

  unsigned short get_number_of_allowed_crossings() const override { return m_number_of_allowed_crossings; }
  void set_number_of_allowed_crossings(unsigned short v) override { m_number_of_allowed_crossings = v; }

  unsigned short get_number_of_tested_crossings() const override { return m_number_of_tested_crossings; }
  void set_number_of_tested_crossings(unsigned short v) override { m_number_of_tested_crossings = v; }

  unsigned short get_number_of_tpc_valid_crossings() const override { return m_number_of_tpc_valid_crossings; }
  void set_number_of_tpc_valid_crossings(unsigned short v) override { m_number_of_tpc_valid_crossings = v; }

  unsigned short get_number_of_vertex_compatible_crossings() const override { return m_number_of_vertex_compatible_crossings; }
  void set_number_of_vertex_compatible_crossings(unsigned short v) override { m_number_of_vertex_compatible_crossings = v; }

  unsigned int get_number_of_candidates() const override { return m_candidates.size(); }
  const TpcCrossingCandidate* get_candidate(unsigned int index) const override { return index < m_candidates.size() ? &m_candidates[index] : nullptr; }
  void add_candidate(const TpcCrossingCandidate& candidate) override { m_candidates.push_back(candidate); }
  void clear_candidates() override { m_candidates.clear(); }

  unsigned char get_status() const override { return m_status; }
  void set_status(unsigned char v) override { m_status = v; }
  void set_status(TpcCrossingStatus status) override { set_status(static_cast<unsigned char>(status)); }

 private:
  static float nan() { return std::numeric_limits<float>::quiet_NaN(); }

  unsigned int m_assembled_track_id {0};
  short m_selected_crossing {std::numeric_limits<short>::max()};
  unsigned int m_silicon_vertex_id {std::numeric_limits<unsigned int>::max()};
  float m_tpc_z0 {nan()};
  float m_silicon_vertex_z {nan()};
  float m_delta_z {nan()};
  float m_best_abs_delta_z {nan()};
  float m_second_best_abs_delta_z {nan()};
  unsigned char m_selected_tier {std::numeric_limits<unsigned char>::max()};
  float m_selected_score {nan()};
  unsigned short m_number_of_available_crossings {0};
  unsigned short m_number_of_allowed_crossings {0};
  unsigned short m_number_of_tested_crossings {0};
  unsigned short m_number_of_tpc_valid_crossings {0};
  unsigned short m_number_of_vertex_compatible_crossings {0};
  std::vector<TpcCrossingCandidate> m_candidates;
  unsigned char m_status {static_cast<unsigned char>(TpcCrossingStatus::Unknown)};

  ClassDefOverride(TpcCrossingDecisionv1, 3)
};

#endif  // TPCTRACKRECO_TPCCROSSINGDECISIONV1_H
