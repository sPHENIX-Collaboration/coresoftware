#include "TowerInfoContainerv2.h"
#include "TowerInfov2.h"

#include <TClonesArray.h>
#include <TSystem.h>

#include <cstdlib>

TowerInfoContainerv2::TowerInfoContainerv2(DETECTOR detec)
  : _detector(detec)
{
  int nchannels = get_channels(detec);
  _clones = new TClonesArray("TowerInfov2", nchannels);
  for (int i = 0; i < nchannels; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerv2::TowerInfoContainerv2(const TowerInfoContainerv2& source)
  : TowerInfoContainer(source)
  , _clones(new TClonesArray("TowerInfov2", source.size()))
  , _detector(source.get_detectorid())
{
  for (unsigned int i = 0; i < source.size(); ++i)
  {
    auto* tower = static_cast<TowerInfov2*>(_clones->ConstructedAt(i, "C"));
    auto* source_tower = static_cast<TowerInfov2*>(source._clones->UncheckedAt(i));
    tower->copy_tower(source_tower);
  }
}

TowerInfoContainerv2::~TowerInfoContainerv2()
{
  delete _clones;
}

void TowerInfoContainerv2::identify(std::ostream& os) const
{
  os << "TowerInfoContainerv2 of size " << size() << std::endl;
}

void TowerInfoContainerv2::Reset()
{
  // clear content of towers in the container for the next event

  for (Int_t i = 0; i < _clones->GetEntriesFast(); ++i)
  {
    TowerInfo* twr = (TowerInfov2*) _clones->UncheckedAt(i);

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

TowerInfov2* TowerInfoContainerv2::get_tower_at_channel(int pos)
{
  return (TowerInfov2*) _clones->At(pos);
}

TowerInfov2* TowerInfoContainerv2::get_tower_at_key(int pos)
{
  int index = decode_key(pos);
  return (TowerInfov2*) _clones->At(index);
}
