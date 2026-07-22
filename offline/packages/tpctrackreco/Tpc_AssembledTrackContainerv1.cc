#include "Tpc_AssembledTrackContainerv1.h"
#include "Tpc_AssembledTrack.h"

ClassImp(Tpc_AssembledTrackContainerv1)

Tpc_AssembledTrackContainerv1::Tpc_AssembledTrackContainerv1()
{
  Reset();
}

Tpc_AssembledTrackContainerv1::~Tpc_AssembledTrackContainerv1()
{
  Reset();
}

void Tpc_AssembledTrackContainerv1::identify(std::ostream& os) const
{
  os << "Tpc_AssembledTrackContainerv1 with "
     << m_tracks.size()
     << " assembled tracks"
     << std::endl;
}

void Tpc_AssembledTrackContainerv1::Reset()
{
  for (auto & m_track : m_tracks)
  {
    delete m_track;
  }
  m_tracks.clear();
}

int Tpc_AssembledTrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* Tpc_AssembledTrackContainerv1::CloneMe() const
{
  Tpc_AssembledTrackContainerv1* copy = new Tpc_AssembledTrackContainerv1();
  for (auto m_track : m_tracks)
  {
    copy->m_tracks.push_back(static_cast<Tpc_AssembledTrack*>(m_track->CloneMe()));
  }
  return copy;
}
