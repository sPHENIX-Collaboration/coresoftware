#ifndef TRACKINFOCONTAINERV1_H
#define TRACKINFOCONTAINERV1_H

#include <phool/PHObject.h>
#include "SvtxTrackInfo.h"
#include "SvtxTrackInfo_v1.h"
#include "TrackInfoContainer.h"

#include <TClonesArray.h>
#include <climits>
#include <cstdint>

class TrackInfoContainer_v1 : public TrackInfoContainer
{
 public:
  TrackInfoContainer_v1();
  ~TrackInfoContainer_v1() override;
  void identify(std::ostream &os = std::cout) const override;
  void Reset() override;

  size_t size() const override { return _clones->GetEntries(); }

  SvtxTrackInfo *get_trackinfo(int pos) override
  {
    return (SvtxTrackInfo *) _clones->At(pos);
  }

  void add_trackinfo(int pos, SvtxTrackInfo trackinfo) override
  {
    new ((*_clones)[pos]) SvtxTrackInfo_v1;
    SvtxTrackInfo_v1 *info = (SvtxTrackInfo_v1 *) _clones->ConstructedAt(pos);
    info->CopyFrom(trackinfo);
  }

  void add_trackinfo(int pos, SvtxTrackInfo *trackinfo) override
  {
    new ((*_clones)[pos]) SvtxTrackInfo_v1;
    SvtxTrackInfo_v1 *info = (SvtxTrackInfo_v1 *) _clones->ConstructedAt(pos);
    info->CopyFrom(trackinfo);
  }

 protected:
  TClonesArray *_clones = nullptr;

 private:
  ClassDefOverride(TrackInfoContainer_v1, 1);
};

#endif
