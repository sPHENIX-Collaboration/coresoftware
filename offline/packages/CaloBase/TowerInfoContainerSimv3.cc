#include "TowerInfoContainerSimv3.h"
#include "TowerInfoSimv3.h"

#include <TClonesArray.h>
#include <TSystem.h>

#include <cassert>

TowerInfoContainerSimv3::TowerInfoContainerSimv3(DETECTOR detec)
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
  _clones = new TClonesArray("TowerInfoSimv3", nchannels);
  for (int i = 0; i < nchannels; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerSimv3::TowerInfoContainerSimv3(const TowerInfoContainerSimv3& source)
  : TowerInfoContainer(source)
  , _clones(new TClonesArray("TowerInfoSimv3", (int) source.size()))
  , _detector(source.get_detectorid())
{
  for (int i = 0; i < (int) source.size(); ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    auto* tower = static_cast<TowerInfoSimv3*>(_clones->ConstructedAt(i, "C"));
    auto* source_tower = static_cast<TowerInfoSimv3*>(source._clones->UncheckedAt(i));
    tower->copy_tower(source_tower);
  }
}

TowerInfoContainerSimv3::~TowerInfoContainerSimv3()
{
  delete _clones;
}

void TowerInfoContainerSimv3::identify(std::ostream& os) const
{
  os << "TowerInfoContainerSimv3 of size " << size() << std::endl;
}

void TowerInfoContainerSimv3::Reset()
{
  // clear content of towers in the container for the next event

  for (Int_t i = 0; i < _clones->GetEntriesFast(); ++i)
  {
    TowerInfo *twr = (TowerInfoSimv3*) _clones->UncheckedAt(i);

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

TowerInfoSimv3* TowerInfoContainerSimv3::get_tower_at_channel(int pos)
{
  return (TowerInfoSimv3*) _clones->At(pos);
}

TowerInfoSimv3* TowerInfoContainerSimv3::get_tower_at_key(int pos)
{
  int index = (int) decode_key(pos);
  return (TowerInfoSimv3*) _clones->At(index);
}
