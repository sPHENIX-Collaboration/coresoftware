#include "Tpc_ModuleTrackContainerv1.h"
#include "Tpc_ModuleTrack.h"

ClassImp(Tpc_ModuleTrackContainerv1)

Tpc_ModuleTrackContainerv1::Tpc_ModuleTrackContainerv1()
{
  Reset();
}

Tpc_ModuleTrackContainerv1::~Tpc_ModuleTrackContainerv1()
{
  Reset();
}

void Tpc_ModuleTrackContainerv1::identify(std::ostream& os) const
{
  os << "Tpc_ModuleTrackContainerv1 with "
     << m_tracks.size()
     << " in-module tracks"
     << std::endl;
}

void Tpc_ModuleTrackContainerv1::Reset()
{
  for (auto* trk : m_tracks)
  {
    delete trk;
  }
  m_tracks.clear();
}

int Tpc_ModuleTrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* Tpc_ModuleTrackContainerv1::CloneMe() const
{
  // Deep copy: clone each owned track.
  auto* copy = new Tpc_ModuleTrackContainerv1();
  for (const auto* trk : m_tracks)
  {
    copy->m_tracks.push_back(static_cast<Tpc_ModuleTrack*>(trk->CloneMe()));
  }
  return copy;
}
