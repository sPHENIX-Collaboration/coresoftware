#include "TrackInfoContainer_v3.h"
#include <phool/PHObject.h>
#include "SvtxTrackInfo.h"
#include "SvtxTrackInfo_v3.h"

#include <climits>
#include <map>

TrackInfoContainer_v3::TrackInfoContainer_v3()
{
  _clones = new TClonesArray("SvtxTrackInfo_v3");
  _clones->SetOwner();
  _clones->SetName("TrackInfoContainer_v3");
}

TrackInfoContainer_v3::~TrackInfoContainer_v3()
{
  _clones->Clear("C");
}

void TrackInfoContainer_v3::identify(std::ostream& os) const
{
  os << "TrackInfoContainer_v3 of size " << size() << std::endl;
}

void TrackInfoContainer_v3::Reset()
{
  _clones->Clear("C");
}
