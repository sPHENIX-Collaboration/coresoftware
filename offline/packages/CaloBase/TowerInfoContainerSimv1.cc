#include "TowerInfoContainerSimv1.h"
#include "TowerInfoSimv1.h"

#include <TClonesArray.h>

#include <cassert>

TowerInfoContainerSimv1::TowerInfoContainerSimv1(DETECTOR detec)
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
    nchannels = 52;
  }
  _clones = new TClonesArray("TowerInfoSimv1", nchannels);
  _clones->SetOwner();
  _clones->SetName("TowerInfoContainerSimv1");
  for (int i = 0; i < nchannels; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerSimv1::TowerInfoContainerSimv1(const TowerInfoContainerSimv1& source)
  : TowerInfoContainer(source)
{
  _detector = source.get_detectorid();
  _clones = new TClonesArray("TowerInfoSimv1", source.size());
  _clones->SetOwner();
  _clones->SetName("TowerInfoContainerSimv1");
  for (unsigned int i = 0; i < source.size(); ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerSimv1::~TowerInfoContainerSimv1()
{
  delete _clones;
}

void TowerInfoContainerSimv1::identify(std::ostream& os) const
{
  os << "TowerInfoContainerSimv1 of size " << size() << std::endl;
}

void TowerInfoContainerSimv1::Reset()
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
}

TowerInfoSimv1* TowerInfoContainerSimv1::get_tower_at_channel(int pos)
{
  return (TowerInfoSimv1*) _clones->At(pos);
}

TowerInfoSimv1* TowerInfoContainerSimv1::get_tower_at_key(int pos)
{
  int index = decode_key(pos);
  return (TowerInfoSimv1*) _clones->At(index);
}

unsigned int TowerInfoContainerSimv1::encode_key(unsigned int towerIndex)
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

unsigned int TowerInfoContainerSimv1::decode_key(unsigned int tower_key)
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
