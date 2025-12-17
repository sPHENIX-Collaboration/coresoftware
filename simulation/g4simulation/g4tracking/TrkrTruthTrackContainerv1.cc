#include "TrkrTruthTrackContainerv1.h"
#include "TrkrTruthTrackv1.h"

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/phool.h>

#include <algorithm>
#include <iostream>
#include <set>

void TrkrTruthTrackContainerv1::Reset()
{
  for (auto entry : m_data)
  {
    delete entry.second;
  }
  m_data.clear();
}

void TrkrTruthTrackContainerv1::addTruthTrack(TrkrTruthTrack* track)
{
  if (!track)
  {
    return;
  }
  const auto id{track->getTrackid()};
  if (hasTrackid(id))
  {
    std::cout << PHWHERE << "Warning, replacing existing track-id(" << id << ")" << std::endl;
  }
  m_data[id] = track;
}

TrkrTruthTrack* TrkrTruthTrackContainerv1::getTruthTrack(unsigned int trackid)
{
  if (!hasTrackid(trackid))
  {
    std::cout << PHWHERE << " Asking for TrkrTruthTrack " << trackid
              << " which is not present. Returning empty track." << std::endl;
    TrkrTruthTrack* track = new TrkrTruthTrackv1();
    return track;
  }
  return m_data[trackid];
}

TrkrTruthTrack* TrkrTruthTrackContainerv1::getTruthTrack(unsigned int id, PHG4TruthInfoContainer* truth_info)
{
  // return the track if already in m_data, otherwise make it and return the newly made track
  if (hasTrackid(id))
  {
    return m_data[id];
  }
  PHG4Particle* particle = /*(PHG4Particlev3*)*/ truth_info->GetParticle(id);
  if (particle == nullptr)
  {
    std::cout << PHWHERE << " Note: embedded track from PHG4TruthInfoContainer, id( "
              << id << " )without an associated PHG4Particle" << std::endl;
    auto *current_track = new TrkrTruthTrackv1();
    current_track->setTrackid(id);
    m_data[id] = current_track;
    return current_track;
  }
  int vtxid = particle->get_vtx_id();
  PHG4VtxPoint* vtx = truth_info->GetVtx(vtxid);
  auto *current_track = new TrkrTruthTrackv1(id, particle, vtx);
  m_data[id] = current_track;
  return current_track;
}

TrkrTruthTrackContainer::ConstRange TrkrTruthTrackContainerv1::getTruthTrackRange() const
{
  return std::make_pair(m_data.begin(), m_data.end());
}

bool TrkrTruthTrackContainerv1::hasTrackid(unsigned int id) const
{
  return (m_data.contains(id));
}

TrkrTruthTrackContainer::Map& TrkrTruthTrackContainerv1::getMap()
{
  return m_data;
}

void TrkrTruthTrackContainerv1::identify(std::ostream& os) const
{
  os << " TrkrTruthTrackContainer data.  Containter " << (int) m_data.size() << " tracks" << std::endl;
  int cnt = 0;
  for (const auto& entry : m_data)
  {
    os << " Track(" << cnt << "): " << std::endl;
    entry.second->identify(os);
    ++cnt;
  }
}
