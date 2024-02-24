#ifndef TRACKINFOCONTAINERV3_H
#define TRACKINFOCONTAINERV3_H

#include "TrackInfoContainer.h"
#include "SvtxTrackInfo_v3.h"
#include "SvtxTrackInfo.h"
#include <phool/PHObject.h>

#include <climits>
#include <map>
#include <TClonesArray.h>
#include <cstdint>

class TrackInfoContainer_v3 : public TrackInfoContainer
{
    public:

    TrackInfoContainer_v3();
    ~TrackInfoContainer_v3() override;
    void identify(std::ostream& os = std::cout) const override;

      void Reset() override;
     //virtual SvtxTrackInfo* get_trackinfo(int /*index*/) override { return nullptr; }
    //virtual TrackInfo* get_tower_at_key(int /*key*/) { return nullptr; }
    //virtual size_t size() const { return 0; }

    size_t size() const override { return _clones->GetEntries(); }

    SvtxTrackInfo_v3* get_trackinfo(int pos) override{
    return (SvtxTrackInfo_v3*) _clones->At(pos);
    }

    void add_trackinfo(int pos, SvtxTrackInfo_v3 trackinfo)
    {
      new((*_clones)[pos]) SvtxTrackInfo_v3;
      SvtxTrackInfo_v3 *info = (SvtxTrackInfo_v3 *)_clones->ConstructedAt(pos);
      info->CopyFrom(trackinfo);

    }

    void add_trackinfo(int pos, SvtxTrackInfo_v3* trackinfo)
    {
      new((*_clones)[pos]) SvtxTrackInfo_v3;
      SvtxTrackInfo_v3 *info = (SvtxTrackInfo_v3 *)_clones->ConstructedAt(pos);
      info->CopyFrom(trackinfo);

    }


    protected:
    TClonesArray *_clones = nullptr;

    private:
    ClassDefOverride(TrackInfoContainer_v3, 1);
};

#endif
