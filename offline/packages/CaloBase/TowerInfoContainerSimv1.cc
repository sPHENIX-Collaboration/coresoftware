#include "TowerInfoContainerSimv1.h"
#include "TowerInfoSimv1.h"

#include <TClonesArray.h>
#include <TSystem.h>

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
  for (int i = 0; i < nchannels; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerSimv1::TowerInfoContainerSimv1(const TowerInfoContainerSimv1& source)
  : TowerInfoContainer(source)
  , _clones(new TClonesArray("TowerInfoSimv1", source.size()))
  , _detector(source.get_detectorid())
{
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
    TowerInfo *twr = (TowerInfoSimv1 *) _clones->UncheckedAt(i);
    
    if (twr == nullptr)
    {
      std::cout << __PRETTY_FUNCTION__ << " Fatal access error:"
                << " _clones->GetSize() = " << _clones->GetSize()
                << " _clones->GetEntriesFast() = " << _clones->GetEntriesFast()
                << " i = " << i << std::endl;
      _clones->Print();
      gSystem->Exit(1);
      exit(1);
    }
    twr->Reset();
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
