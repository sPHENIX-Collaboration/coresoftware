#include "Tpc_PolyTrackVertexv1.h"

#include <cmath>

ClassImp(Tpc_PolyTrackVertexv1)

Tpc_PolyTrackVertexv1::Tpc_PolyTrackVertexv1()
{
  Reset();
}

void Tpc_PolyTrackVertexv1::identify(std::ostream& os) const
{
  os << "Tpc_PolyTrackVertexv1:"
     << " track_id=" << m_track_id
     << " source_assembled_track_id=" << m_source_assembled_track_id
     << " dca2d=" << m_dca2d
     << " z0=" << m_z0
     << " pca_valid=" << m_pca_valid
     << " pca=(" << m_pca_x << ", " << m_pca_y << ", " << m_pca_z << ")"
     << std::endl;
}

void Tpc_PolyTrackVertexv1::Reset()
{
  m_track_id = 0;
  m_source_assembled_track_id = 0;
  m_dca2d = 0.0;
  m_z0 = 0.0;
  m_pca_valid = 0;
  m_pca_x = 0.0;
  m_pca_y = 0.0;
  m_pca_z = 0.0;
  m_pca_radius = 0.0;
  m_pca_phi = 0.0;
}

int Tpc_PolyTrackVertexv1::isValid() const
{
  return m_pca_valid && std::isfinite(m_z0);
}
