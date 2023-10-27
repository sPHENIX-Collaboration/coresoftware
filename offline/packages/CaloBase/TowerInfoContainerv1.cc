#include "TowerInfoContainerv1.h"
#include "TowerInfov1.h"
#include "TowerInfoDefs.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <TClonesArray.h>

#include <cassert>


TowerInfoContainerv1::TowerInfoContainerv1(DETECTOR detec)
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
  _clones = new TClonesArray("TowerInfov1", nchannels);
  _clones->SetOwner();
  _clones->SetName("TowerInfoContainerv1");
  for (int i = 0; i < nchannels; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerv1::~TowerInfoContainerv1()
{
  delete _clones;
}

void TowerInfoContainerv1::identify(std::ostream& os) const
{
  os << "TowerInfoContainerv1 of size " << size() << std::endl;
}

void TowerInfoContainerv1::Reset()
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

TowerInfov1* TowerInfoContainerv1::get_tower_at_channel(int pos)
{
  return (TowerInfov1*) _clones->At(pos);
}


TowerInfov1* TowerInfoContainerv1::get_tower_at_key(int pos)
{
  int index = decode_key(pos);
  return (TowerInfov1*) _clones->At(index);
}

unsigned int TowerInfoContainerv1::encode_key(unsigned int towerIndex)
{
  int key = 0;
  if (_detector == DETECTOR::EMCAL)
    {
      key = TowerInfoContainerv1::encode_emcal(towerIndex);
    }
  else if (_detector == DETECTOR::HCAL)
    {
      key = TowerInfoContainerv1::encode_hcal(towerIndex);
    }
  else if (_detector == DETECTOR::SEPD)
    {
    key = TowerInfoContainerv1::encode_epd(towerIndex);
    }
  else if (_detector == DETECTOR::MBD)
    {
    key = TowerInfoContainerv1::encode_mbd(towerIndex);
    }
  else if (_detector == DETECTOR::ZDC)
    {
    key = TowerInfoContainerv1::encode_zdc(towerIndex);
    }
  return key;
}

unsigned int TowerInfoContainerv1::decode_key(unsigned int tower_key)
{
  int index = 0;

  if (_detector == DETECTOR::EMCAL)
  {
    index = TowerInfoContainerv1::decode_emcal(tower_key);
  }
  else if (_detector == DETECTOR::HCAL)
  {
    index = TowerInfoContainerv1::decode_hcal(tower_key);
  }
  else if (_detector == DETECTOR::SEPD)
  {
    index = TowerInfoContainerv1::decode_epd(tower_key);
  }
  else if (_detector == DETECTOR::MBD)
  {
    index = TowerInfoContainerv1::decode_mbd(tower_key);
  }
  else if (_detector == DETECTOR::ZDC)
  {
    index = TowerInfoContainerv1::decode_zdc(tower_key);
  }
  return index;
}


