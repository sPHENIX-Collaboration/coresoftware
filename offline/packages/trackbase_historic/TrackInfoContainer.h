#ifndef TRACKINFOCONTAINER_H
#define TRACKINFOCONTAINER_H

#include "SvtxTrackInfo.h"
#include "SvtxTrackInfo_v1.h"

#include <phool/PHObject.h>

#include <climits>
#include <map>

class TrackInfoContainer : public PHObject
{
 public:
  TrackInfoContainer() = default;
  ~TrackInfoContainer() override = default;
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrackInfoContainer base class" << std::endl;
  }

  virtual void Reset() override {}
  virtual SvtxTrackInfo* get_trackinfo(int /*index*/) { return nullptr; }
  virtual void add_trackinfo(int, SvtxTrackInfo) {}
  virtual void add_trackinfo(int, SvtxTrackInfo*) {}
  // virtual TrackInfo* get_tower_at_key(int /*key*/) { return nullptr; }
  virtual size_t size() const { return 0; }

 private:
  ClassDefOverride(TrackInfoContainer, 1);
};

#endif