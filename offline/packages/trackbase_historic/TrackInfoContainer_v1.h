#ifndef TRACKINFOCONTAINERV1_H
#define TRACKINFOCONTAINERV1_H

#include "TrackInfoContainer.h"
#include "SvtxTrackInfo_v1.h"
#include "SvtxTrackInfo.h"
#include <phool/PHObject.h>

#include <climits>
#include <map>
#include <TClonesArray.h>
#include <cstdint>

class TrackInfoContainer_v1 : public TrackInfoContainer
{
    public:
 
    TrackInfoContainer_v1();
    ~TrackInfoContainer_v1() override;
    void identify(std::ostream& os = std::cout) const override;

      void Reset() override;
     //virtual SvtxTrackInfo* get_trackinfo(int /*index*/) override { return nullptr; }
    //virtual TrackInfo* get_tower_at_key(int /*key*/) { return nullptr; }
    //virtual size_t size() const { return 0; }

    size_t size() const override { return _clones->GetEntries(); }

    SvtxTrackInfo* get_trackinfo(int pos) override{
    return (SvtxTrackInfo*) _clones->At(pos);
    }

    //    void TrackInfoContainer_v1::add_trackinfo(int, SvtxTrackInfo_v1 ) {}
    void add_trackinfo(int pos, SvtxTrackInfo trackinfo) override{
      new((*_clones)[pos]) SvtxTrackInfo;
      SvtxTrackInfo *info = (SvtxTrackInfo *)_clones->ConstructedAt(pos);
      info->CopyFrom(trackinfo);

    }

    void add_trackinfo(int pos, SvtxTrackInfo* trackinfo) override{
      new((*_clones)[pos]) SvtxTrackInfo;
      SvtxTrackInfo *info = (SvtxTrackInfo *)_clones->ConstructedAt(pos);
      info->CopyFrom(trackinfo);

    }


    protected:
    TClonesArray *_clones = nullptr;

    private:
    ClassDefOverride(TrackInfoContainer_v1, 1);
};

#endif