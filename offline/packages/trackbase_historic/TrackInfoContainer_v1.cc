#include "SvtxTrackInfo.h"
#include "TrackInfoContainer_v1.h"
#include <phool/PHObject.h>

#include <climits>
#include <map>


TrackInfoContainer_v1::TrackInfoContainer_v1()
{
  _clones = new TClonesArray("TrackInfoContainer_v1");
  _clones->SetOwner();
  _clones->SetName("TrackInfoContainer_v1");
}

TrackInfoContainer_v1::~TrackInfoContainer_v1()
{
  _clones->Clear("C");
}

void TrackInfoContainer_v1::identify(std::ostream& os) const
{
  os << "TrackInfoContainer_v1 of size " << size() << std::endl;
}

void TrackInfoContainer_v1::Reset()
{
  _clones->Clear("C");
}
