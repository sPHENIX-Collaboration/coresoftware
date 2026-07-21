#pragma once

#include "Tpc_AssembledTrackContainer.h"

#include <iostream>
#include <vector>

class Tpc_AssembledTrack;

class Tpc_AssembledTrackContainerv1 : public Tpc_AssembledTrackContainer
{
 public:
  Tpc_AssembledTrackContainerv1();
  ~Tpc_AssembledTrackContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override
  {
    return static_cast<unsigned int>(m_tracks.size());
  }

  void add_track(Tpc_AssembledTrack* trk) override
  {
    m_tracks.push_back(trk);
  }

  const Tpc_AssembledTrack* get_track(unsigned int i) const override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

  Tpc_AssembledTrack* get_track(unsigned int i) override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

 private:
  std::vector<Tpc_AssembledTrack*> m_tracks;

  ClassDefOverride(Tpc_AssembledTrackContainerv1, 1)
};
