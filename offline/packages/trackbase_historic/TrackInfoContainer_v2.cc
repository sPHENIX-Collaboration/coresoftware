#include "TrackInfoContainer_v2.h"
#include <phool/PHObject.h>
#include "SvtxTrackInfo_v2.h"

TrackInfoContainer_v2::TrackInfoContainer_v2()
{
  _clones = new TClonesArray("SvtxTrackInfo_v2");
  _clones->SetOwner();
  _clones->SetName("TrackInfoContainer_v2");
}

TrackInfoContainer_v2::~TrackInfoContainer_v2()
{
  _clones->Clear("C");
}

void TrackInfoContainer_v2::identify(std::ostream& os) const
{
  os << "TrackInfoContainer_v2 of size " << size() << std::endl;
}

void TrackInfoContainer_v2::Reset()
{
  _clones->Clear("C");
}
