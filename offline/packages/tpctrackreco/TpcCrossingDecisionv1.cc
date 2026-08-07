#include "TpcCrossingDecisionv1.h"

ClassImp(TpcCrossingCandidate)
ClassImp(TpcCrossingDecisionv1)

TpcCrossingDecisionv1::TpcCrossingDecisionv1()
{
  Reset();
}

void TpcCrossingDecisionv1::identify(std::ostream& os) const
{
  os << "TpcCrossingDecisionv1:"
     << " assembled_track_id=" << m_assembled_track_id
     << " selected_crossing=" << m_selected_crossing
     << " status=" << static_cast<unsigned int>(m_status)
     << " selected_tier=" << static_cast<unsigned int>(m_selected_tier)
     << " selected_score=" << m_selected_score
     << " tpc_z0=" << m_tpc_z0
     << " silicon_vertex_id=" << m_silicon_vertex_id
     << " delta_z=" << m_delta_z
     << " available=" << m_number_of_available_crossings
     << " allowed=" << m_number_of_allowed_crossings
     << " tested=" << m_number_of_tested_crossings
     << " tpc_valid=" << m_number_of_tpc_valid_crossings
     << " vertex_compatible=" << m_number_of_vertex_compatible_crossings
     << " candidates=" << m_candidates.size()
     << std::endl;
}

void TpcCrossingDecisionv1::Reset()
{
  m_assembled_track_id = 0;
  m_selected_crossing = std::numeric_limits<short>::max();
  m_silicon_vertex_id = std::numeric_limits<unsigned int>::max();
  m_tpc_z0 = nan();
  m_silicon_vertex_z = nan();
  m_delta_z = nan();
  m_best_abs_delta_z = nan();
  m_second_best_abs_delta_z = nan();
  m_selected_tier = std::numeric_limits<unsigned char>::max();
  m_selected_score = nan();
  m_number_of_available_crossings = 0;
  m_number_of_allowed_crossings = 0;
  m_number_of_tested_crossings = 0;
  m_number_of_tpc_valid_crossings = 0;
  m_number_of_vertex_compatible_crossings = 0;
  m_candidates.clear();
  m_status = static_cast<unsigned char>(TpcCrossingStatus::Unknown);
}

int TpcCrossingDecisionv1::isValid() const
{
  return m_status != static_cast<unsigned char>(TpcCrossingStatus::Unknown) ? 1 : 0;
}
