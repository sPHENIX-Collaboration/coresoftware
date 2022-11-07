#include "TrkrTruthTrackContainerv1.h"
#include "TrkrTruthTrackv1.h"

#include <algorithm>

void TrkrTruthTrackContainerv1::Reset()
{
  for (auto& track : m_data)
  {
    delete track;
  }
  m_data.clear();
}

void TrkrTruthTrackContainerv1::addTruthTrack(TrkrTruthTrack* track)
{
  // note: assumed that track is ordered numerically by trackid already
  m_data.push_back(track);
}

TrkrTruthTrackContainer::ConstRange TrkrTruthTrackContainerv1::getTruthTrackRange() const
{
  return std::make_pair(m_data.begin(), m_data.end());
}

bool TrkrTruthTrackContainerv1::hasTrackid(unsigned int trackid) const
{
  return std::binary_search(m_data.begin(), m_data.end(), trackid, TrkrTruthTrack::Comp());
}

TrkrTruthTrackContainerv1::Vector& TrkrTruthTrackContainerv1::getTruthTracks()
{
  return m_data;
}

TrkrTruthTrack* TrkrTruthTrackContainerv1::getTruthTrack(unsigned int trackid) const
{
  auto iter = std::lower_bound(m_data.begin(), m_data.end(), trackid, TrkrTruthTrack::Comp());
  if (iter == m_data.end())
  {
    std::cout << "Asking for TrkrTruthTrack " << trackid << " which is not present. Returning empty track." << std::endl;
    TrkrTruthTrack* track = new TrkrTruthTrackv1();
    return track;
  }
  return *iter;
}

void TrkrTruthTrackContainerv1::identify(std::ostream& os) const
{
  os << " TrkrTruthTrackContainer data.  Containter " << m_data.size() << " tracks" << std::endl;
  int cnt = 0;
  for (auto& entry : m_data)
  {
    os << " Track(" << cnt << "): " << std::endl;
    entry->identify(os);
    ++cnt;
  }
}
