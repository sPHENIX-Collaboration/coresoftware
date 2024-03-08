#ifndef TRACKINFOCONTAINERV2_H
#define TRACKINFOCONTAINERV2_H

#include <phool/PHObject.h>
#include "SvtxTrackInfo_v2.h"
#include "TrackInfoContainer.h"

#include <TClonesArray.h>

class TrackInfoContainer_v2 : public TrackInfoContainer
{
 public:
  TrackInfoContainer_v2();
  ~TrackInfoContainer_v2() override;
  void identify(std::ostream &os = std::cout) const override;
  void Reset() override;

  size_t size() const override { return _clones->GetEntries(); }

  SvtxTrackInfo_v2 *get_trackinfo(int pos) override
  {
    return (SvtxTrackInfo_v2 *) _clones->At(pos);
  }

  void add_trackinfo(int pos, SvtxTrackInfo trackinfo) override
  {
    new ((*_clones)[pos]) SvtxTrackInfo_v2;
    SvtxTrackInfo_v2 *info = (SvtxTrackInfo_v2 *) _clones->ConstructedAt(pos);
    info->CopyFrom(trackinfo);
  }

  void add_trackinfo(int pos, SvtxTrackInfo *trackinfo) override
  {
    new ((*_clones)[pos]) SvtxTrackInfo_v2;
    SvtxTrackInfo_v2 *info = (SvtxTrackInfo_v2 *) _clones->ConstructedAt(pos);
    info->CopyFrom(trackinfo);
  }

 protected:
  TClonesArray *_clones = nullptr;

 private:
  ClassDefOverride(TrackInfoContainer_v2, 1);
};

#endif
