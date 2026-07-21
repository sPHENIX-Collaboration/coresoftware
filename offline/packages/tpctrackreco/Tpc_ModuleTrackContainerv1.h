#pragma once

#include "Tpc_ModuleTrackContainer.h"

#include <iostream>
#include <vector>

class Tpc_ModuleTrack;

class Tpc_ModuleTrackContainerv1 : public Tpc_ModuleTrackContainer
{
 public:
  Tpc_ModuleTrackContainerv1();
  ~Tpc_ModuleTrackContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override
  {
    return static_cast<unsigned int>(m_tracks.size());
  }

  // Takes ownership of the track pointer.
  void add_track(Tpc_ModuleTrack* trk) override
  {
    m_tracks.push_back(trk);
  }

  const Tpc_ModuleTrack* get_track(unsigned int i) const override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

  Tpc_ModuleTrack* get_track(unsigned int i) override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

 private:
  std::vector<Tpc_ModuleTrack*> m_tracks;

  ClassDefOverride(Tpc_ModuleTrackContainerv1, 1)
};
