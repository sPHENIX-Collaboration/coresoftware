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
  for (auto& m_track : m_tracks)
  {
    delete m_track;
  }
  m_tracks.clear();
}

int Tpc_PolyTrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* Tpc_PolyTrackContainerv1::CloneMe() const
{
  Tpc_PolyTrackContainerv1* copy = new Tpc_PolyTrackContainerv1();
  for (auto* m_track : m_tracks)
  {
    if (m_track)
    {
      copy->m_tracks.push_back(static_cast<Tpc_PolyTrack*>(m_track->CloneMe()));
    }
  }
  return copy;
}
