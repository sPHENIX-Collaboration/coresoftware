#pragma once

#include "Tpc_PolyTrackContainer.h"

#include <iostream>
#include <vector>

class Tpc_PolyTrack;

class Tpc_PolyTrackContainerv1 : public Tpc_PolyTrackContainer
{
 public:
  Tpc_PolyTrackContainerv1();
  ~Tpc_PolyTrackContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_tracks.size()); }
  void add_track(Tpc_PolyTrack* trk) override { m_tracks.push_back(trk); }
  const Tpc_PolyTrack* get_track(unsigned int i) const override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }
  Tpc_PolyTrack* get_track(unsigned int i) override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

 private:
  std::vector<Tpc_PolyTrack*> m_tracks;

  ClassDefOverride(Tpc_PolyTrackContainerv1, 1)
};
