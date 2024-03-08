#ifndef TRACKINFOCONTAINERV3_H
#define TRACKINFOCONTAINERV3_H

#include <phool/PHObject.h>
#include "SvtxTrackInfo.h"
#include "SvtxTrackInfo_v3.h"
#include "TrackInfoContainer.h"

#include <TClonesArray.h>

class TrackInfoContainer_v3 : public TrackInfoContainer
{
 public:
  TrackInfoContainer_v3();
  ~TrackInfoContainer_v3() override;
  void identify(std::ostream &os = std::cout) const override;
  void Reset() override;

  size_t size() const override { return _clones->GetEntries(); }

  SvtxTrackInfo_v3 *get_trackinfo(int pos) override
  {
    return (SvtxTrackInfo_v3 *) _clones->At(pos);
  }

  void add_trackinfo(int pos, SvtxTrackInfo trackinfo) override
  {
    new ((*_clones)[pos]) SvtxTrackInfo_v3;
    SvtxTrackInfo_v3 *info = (SvtxTrackInfo_v3 *) _clones->ConstructedAt(pos);
    info->CopyFrom(trackinfo);
  }

  void add_trackinfo(int pos, SvtxTrackInfo *trackinfo) override
  {
    new ((*_clones)[pos]) SvtxTrackInfo_v3;
    SvtxTrackInfo_v3 *info = (SvtxTrackInfo_v3 *) _clones->ConstructedAt(pos);
    info->CopyFrom(trackinfo);
  }

 protected:
  TClonesArray *_clones = nullptr;

 private:
  ClassDefOverride(TrackInfoContainer_v3, 1);
};

#endif
