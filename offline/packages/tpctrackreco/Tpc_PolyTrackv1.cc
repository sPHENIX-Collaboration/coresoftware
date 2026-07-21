#include "Tpc_PolyTrackv1.h"

#include <cmath>
#include <limits>

ClassImp(Tpc_PolyTrackv1)

Tpc_PolyTrackv1::Tpc_PolyTrackv1()
{
  Reset();
}

void Tpc_PolyTrackv1::identify(std::ostream& os) const
{
  os << "Tpc_PolyTrackv1:"
     << " event=" << m_event
     << " track_id=" << m_track_id
     << " source_assembled_track_id=" << m_source_assembled_track_id
     << " fit_status=" << m_fit_status
     << " nclusters=" << m_nclusters
     << " pos=(" << m_x << "," << m_y << "," << m_z << ")"
     << " mom=(" << m_px << "," << m_py << "," << m_pz << ")"
     << " charge=" << m_charge
     << " chi2=" << m_chi2
     << " ndf=" << m_ndf
     << " dedx=" << m_dedx
     << std::endl;
}

void Tpc_PolyTrackv1::Reset()
{
  m_event = 0;
  m_track_id = 0;
  m_source_assembled_track_id = 0;
  m_fit_status = 0;
  m_nclusters = 0;
  m_x = 0.0;
  m_y = 0.0;
  m_z = 0.0;
  m_px = 0.0;
  m_py = 0.0;
  m_pz = 0.0;
  m_charge = 0.0;
  m_chi2 = 0.0;
  m_ndf = 0.0;
  m_dedx = std::numeric_limits<double>::quiet_NaN();
  m_cov.assign(36, 0.0);
}

int Tpc_PolyTrackv1::isValid() const
{
  return (m_fit_status != 0 && m_nclusters > 0 &&
          std::isfinite(m_x) && std::isfinite(m_y) && std::isfinite(m_z) &&
          std::isfinite(m_px) && std::isfinite(m_py) && std::isfinite(m_pz)) ? 1 : 0;
}
