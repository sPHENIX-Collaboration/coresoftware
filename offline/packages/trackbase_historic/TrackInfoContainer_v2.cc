#include "SvtxTrackInfo_v2.h"
#include "TrackInfoContainer_v2.h"
#include <phool/PHObject.h>

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
  // clear content of towers in the container for the next event
  for (Int_t i = 0; i < _clones->GetEntriesFast(); ++i)
  {
    TObject* obj = _clones->UncheckedAt(i);

    if (obj == nullptr)
    {
      std::cout << __PRETTY_FUNCTION__ << " Fatal access error:"
                << " _clones->GetSize() = " << _clones->GetSize()
                << " _clones->GetEntriesFast() = " << _clones->GetEntriesFast()
                << " i = " << i << std::endl;
      _clones->Print();
    }

    assert(obj);
    // same as TClonesArray::Clear() but only clear but not to erase all towers
    obj->Clear();
    obj->ResetBit(kHasUUID);
    obj->ResetBit(kIsReferenced);
    obj->SetUniqueID(0);
  }
  _clones->Clear("C");
}
