#include "Tpc_ModuleTrackv1.h"

ClassImp(Tpc_ModuleTrackv1)

Tpc_ModuleTrackv1::Tpc_ModuleTrackv1()
{
  Reset();
}

void Tpc_ModuleTrackv1::identify(std::ostream& os) const
{
  os << "Tpc_ModuleTrackv1:"
     << " event=" << m_event
     << " track_id=" << m_track_id
     << " region=" << m_region
     << " sector=" << m_sector
     << " side=" << m_side
     << " nblobs=" << m_nblobs
     << " nrawhits=" << m_nrawhits
     << " first_layer=" << m_first_layer
     << " last_layer=" << m_last_layer
     << " nhit_indices=" << m_hit_indices.size()
     << std::endl;
}

void Tpc_ModuleTrackv1::Reset()
{
  m_event = 0;
  m_track_id = 0;

  m_region = 0;
  m_sector = 0;
  m_side = 0;

  m_nblobs = 0;
  m_nrawhits = 0;

  m_first_layer = 0;
  m_last_layer = 0;

  m_hit_indices.clear();
}

int Tpc_ModuleTrackv1::isValid() const
{
  return m_hit_indices.empty() ? 0 : 1;
}
