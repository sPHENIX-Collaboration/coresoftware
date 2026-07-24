#pragma once

#include <phool/PHObject.h>

#include <iostream>

class Tpc_AssembledTrack;

class Tpc_AssembledTrackContainer : public PHObject
{
 public:
  Tpc_AssembledTrackContainer() = default;
  ~Tpc_AssembledTrackContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_AssembledTrackContainer base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  virtual unsigned int size() const { return 0; }
  virtual void add_track(Tpc_AssembledTrack*) {}
  virtual const Tpc_AssembledTrack* get_track(unsigned int) const { return nullptr; }
  virtual Tpc_AssembledTrack* get_track(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(Tpc_AssembledTrackContainer, 0)
};
