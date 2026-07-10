#include "TowerInfoContainerv3.h"
#include "TowerInfov3.h"

#include <TClonesArray.h>
#include <TSystem.h>

#include <cstdlib>

TowerInfoContainerv3::TowerInfoContainerv3(DETECTOR detec)
  : _detector(detec)
{
  int nchannels = get_channels(detec);
  _clones = new TClonesArray("TowerInfov3", nchannels);
  for (int i = 0; i < nchannels; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerv3::TowerInfoContainerv3(const TowerInfoContainerv3& source)
  : TowerInfoContainer(source)
  , _clones(new TClonesArray("TowerInfov3", source.size()))
  , _detector(source.get_detectorid())
{
  for (unsigned int i = 0; i < source.size(); ++i)
  {
    auto* tower = static_cast<TowerInfov3*>(_clones->ConstructedAt(i, "C"));
    auto* source_tower = static_cast<TowerInfov3*>(source._clones->UncheckedAt(i));
    tower->copy_tower(source_tower);
  }
}

TowerInfoContainerv3::~TowerInfoContainerv3()
{
  delete _clones;
}

void TowerInfoContainerv3::identify(std::ostream& os) const
{
  os << "TowerInfoContainerv3 of size " << size() << std::endl;
}

void TowerInfoContainerv3::Reset()
{
  // clear content of towers in the container for the next event

  for (Int_t i = 0; i < _clones->GetEntriesFast(); ++i)
  {
    TowerInfo* twr = (TowerInfov3*) _clones->UncheckedAt(i);

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

TowerInfov3* TowerInfoContainerv3::get_tower_at_channel(int pos)
{
  return (TowerInfov3*) _clones->At(pos);
}

TowerInfov3* TowerInfoContainerv3::get_tower_at_key(int pos)
{
  int index = decode_key(pos);
  return (TowerInfov3*) _clones->At(index);
}
