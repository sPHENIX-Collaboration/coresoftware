#pragma once

#include <phool/PHObject.h>

#include <iostream>

class Tpc_PolyTrack;

class Tpc_PolyTrackContainer : public PHObject
{
 public:
  Tpc_PolyTrackContainer() = default;
  ~Tpc_PolyTrackContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_PolyTrackContainer base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  virtual unsigned int size() const { return 0; }
  virtual void add_track(Tpc_PolyTrack*) {}
  virtual const Tpc_PolyTrack* get_track(unsigned int) const { return nullptr; }
  virtual Tpc_PolyTrack* get_track(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(Tpc_PolyTrackContainer, 0)
};
