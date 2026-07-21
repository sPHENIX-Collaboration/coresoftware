#pragma once

#include <phool/PHObject.h>

#include <iostream>

class Tpc_ModuleTrack;

// Abstract container for in-module tracks.
// Concrete storage is in versioned subclasses (Tpc_ModuleTrackContainerv1, ...).
class Tpc_ModuleTrackContainer : public PHObject
{
 public:
  Tpc_ModuleTrackContainer() = default;
  ~Tpc_ModuleTrackContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_ModuleTrackContainer base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  virtual unsigned int size() const { return 0; }

  virtual void add_track(Tpc_ModuleTrack* /*trk*/) {}

  virtual const Tpc_ModuleTrack* get_track(unsigned int /*i*/) const { return nullptr; }
  virtual Tpc_ModuleTrack* get_track(unsigned int /*i*/) { return nullptr; }

 private:
  ClassDefOverride(Tpc_ModuleTrackContainer, 0)
};
