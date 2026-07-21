#include "Tpc_PolyClusterv1.h"

#include <cmath>

ClassImp(Tpc_PolyClusterv1)

Tpc_PolyClusterv1::Tpc_PolyClusterv1()
{
  Reset();
}

void Tpc_PolyClusterv1::identify(std::ostream& os) const
{
  os << "Tpc_PolyClusterv1:"
     << " event=" << m_event
     << " cluster_id=" << m_cluster_id
     << " source_assembled_track_id=" << m_source_assembled_track_id
     << " side=" << m_side
     << " nhits=" << m_hit_indices.size()
     << " centroid_x=" << m_centroid_x
     << " centroid_y=" << m_centroid_y
     << " centroid_z=" << m_centroid_z
     << " rms_x=" << m_rms_x
     << " rms_y=" << m_rms_y
     << " rms_z=" << m_rms_z
     << " adc=" << m_adc
     << " phi_width=" << m_phi_width
     << " time_width=" << m_time_width
     << " phase=" << m_phase
     << std::endl;
}

void Tpc_PolyClusterv1::Reset()
{
  m_event = 0;
  m_cluster_id = 0;
  m_source_assembled_track_id = 0;
  m_side = 0;
  m_centroid_x = 0.0;
  m_centroid_y = 0.0;
  m_centroid_z = 0.0;
  m_rms_x = 0.0;
  m_rms_y = 0.0;
  m_rms_z = 0.0;
  m_adc = 0.0;
  m_phi_width = 0;
  m_time_width = 0;
  m_phase = 0.0;
  m_hit_indices.clear();
  m_hit_x.clear();
  m_hit_y.clear();
  m_hit_z.clear();
}

int Tpc_PolyClusterv1::isValid() const
{
  return (!m_hit_indices.empty() &&
          std::isfinite(m_centroid_x) &&
          std::isfinite(m_centroid_y) &&
          std::isfinite(m_centroid_z)) ? 1 : 0;
}
