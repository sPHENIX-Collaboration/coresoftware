#include "BbcPmtInfoContainerV1.h"
#include "BbcPmtInfoV1.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <TClonesArray.h>

#include <cassert>


BbcPmtInfoContainerV1::BbcPmtInfoContainerV1()
  : _detector( DETECTOR::MBD )
{
  const int nchannels = 128;
  _clones = new TClonesArray("BbcPmtInfoV1", nchannels);
  _clones->SetOwner();
  _clones->SetName("BbcPmtInfoContainerV1");
  for (int i = 0; i < nchannels; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

BbcPmtInfoContainerV1::~BbcPmtInfoContainerV1()
{
  delete _clones;
}

void BbcPmtInfoContainerV1::Reset()
{
  // clear content of towers in the container for the next event

  for (Int_t i = 0; i < _clones->GetEntriesFast(); ++i)
  {
    TObject* obj = _clones->UncheckedAt(i);

    if (obj==nullptr)
    {
      std::cout<<__PRETTY_FUNCTION__<<" Fatal access error:"
          <<" _clones->GetSize() = "<<_clones->GetSize()
          <<" _clones->GetEntriesFast() = "<<_clones->GetEntriesFast()
          <<" i = "<<i<<std::endl;
      _clones->Print();
    }

    assert(obj);
    // same as TClonesArray::Clear() but only clear but not to erase all towers
    obj->Clear();
    obj->ResetBit(kHasUUID);
    obj->ResetBit(kIsReferenced);
    obj->SetUniqueID(0);
  }
}

BbcPmtInfoV1* BbcPmtInfoContainerV1::get_tower_at_channel(int pos)
{
  return (BbcPmtInfoV1*)_clones->At(pos);
}

