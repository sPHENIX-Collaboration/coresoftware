#include "TowerInfoContainerv2.h"
#include "TowerInfov2.h"
#include "TowerInfoDefs.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <TClonesArray.h>

#include <cassert>


TowerInfoContainerv2::TowerInfoContainerv2(DETECTOR detec)
  : _detector(detec)
{
  int nchannels = 744;
  if (_detector == DETECTOR::SEPD)
  {
    nchannels = 744;
  }
  else if (_detector == DETECTOR::EMCAL)
  {
    nchannels = 24576;
  }
  else if (_detector == DETECTOR::HCAL)
  {
    nchannels = 1536;
  }
  else if (_detector == DETECTOR::MBD)
  {
    nchannels = 256;
  }
  else if (_detector == DETECTOR::ZDC)
  {
    nchannels = 16;
  }
  _clones = new TClonesArray("TowerInfov2", nchannels);
  _clones->SetOwner();
  _clones->SetName("TowerInfoContainerv2");
  for (int i = 0; i < nchannels; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerv2::~TowerInfoContainerv2()
{
  delete _clones;
}

void TowerInfoContainerv2::Reset()
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

TowerInfov2* TowerInfoContainerv2::get_tower_at_channel(int pos)
{
  return (TowerInfov2*) _clones->At(pos);
}


TowerInfov2* TowerInfoContainerv2::get_tower_at_key(int pos)
{
  int index = decode_key(pos);
  return (TowerInfov2*) _clones->At(index);
}

unsigned int TowerInfoContainerv2::encode_key(unsigned int towerIndex)
{
  int key = 0;
  if (_detector == DETECTOR::EMCAL)
    {
      key = TowerInfoContainer::encode_emcal(towerIndex);
    }
  else if (_detector == DETECTOR::HCAL)
    {
      key = TowerInfoContainer::encode_hcal(towerIndex);
    }
  else if (_detector == DETECTOR::SEPD)
    {
    key = TowerInfoContainer::encode_epd(towerIndex);
    }
  else if (_detector == DETECTOR::MBD)
    {
    key = TowerInfoContainer::encode_mbd(towerIndex);
    }
  else if (_detector == DETECTOR::ZDC)
    {
    key = TowerInfoContainer::encode_zdc(towerIndex);
    }
  return key;
}

unsigned int TowerInfoContainerv2::decode_key(unsigned int tower_key)
{
  int index = 0;

  if (_detector == DETECTOR::EMCAL)
  {
    index = TowerInfoContainer::decode_emcal(tower_key);
  }
  else if (_detector == DETECTOR::HCAL)
  {
    index = TowerInfoContainer::decode_hcal(tower_key);
  }
  else if (_detector == DETECTOR::SEPD)
  {
    index = TowerInfoContainer::decode_epd(tower_key);
  }
  else if (_detector == DETECTOR::MBD)
  {
    index = TowerInfoContainer::decode_mbd(tower_key);
  }
  else if (_detector == DETECTOR::ZDC)
  {
    index = TowerInfoContainer::decode_zdc(tower_key);
  }
  return index;
}


