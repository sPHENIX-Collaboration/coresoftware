#include "Tpc_PolyTrackContainerv1.h"
#include "Tpc_PolyTrack.h"

ClassImp(Tpc_PolyTrackContainerv1)

Tpc_PolyTrackContainerv1::Tpc_PolyTrackContainerv1()
{
  Reset();
}

Tpc_PolyTrackContainerv1::~Tpc_PolyTrackContainerv1()
{
  Reset();
}

void Tpc_PolyTrackContainerv1::identify(std::ostream& os) const
{
  os << "Tpc_PolyTrackContainerv1 with " << m_tracks.size() << " tracks" << std::endl;
}

void Tpc_PolyTrackContainerv1::Reset()
{
  for (unsigned int i = 0; i < m_tracks.size(); ++i) delete m_tracks[i];
  m_tracks.clear();
}

int Tpc_PolyTrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* Tpc_PolyTrackContainerv1::CloneMe() const
{
  Tpc_PolyTrackContainerv1* copy = new Tpc_PolyTrackContainerv1();
  for (unsigned int i = 0; i < m_tracks.size(); ++i)
  {
    if (m_tracks[i]) copy->m_tracks.push_back(static_cast<Tpc_PolyTrack*>(m_tracks[i]->CloneMe()));
  }
  return copy;
}
