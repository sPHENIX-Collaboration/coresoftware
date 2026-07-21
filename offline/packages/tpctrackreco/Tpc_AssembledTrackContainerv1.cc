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
  for (unsigned int i = 0; i < m_tracks.size(); ++i)
  {
    delete m_tracks[i];
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
  for (unsigned int i = 0; i < m_tracks.size(); ++i)
  {
    copy->m_tracks.push_back(static_cast<Tpc_AssembledTrack*>(m_tracks[i]->CloneMe()));
  }
  return copy;
}
