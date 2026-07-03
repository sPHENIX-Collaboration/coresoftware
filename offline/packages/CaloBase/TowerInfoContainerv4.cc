#include "TowerInfoContainerv4.h"
#include "TowerInfov4.h"

#include <TClonesArray.h>
#include <TSystem.h>

#include <cassert>

TowerInfoContainerv4::TowerInfoContainerv4(DETECTOR detec)
  : _detector(detec)
{
  int nchannels = get_channels(detec);
  _clones = new TClonesArray("TowerInfov4", nchannels);
  for (int i = 0; i < nchannels; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerv4::TowerInfoContainerv4(const TowerInfoContainerv4& source)
  : TowerInfoContainer(source)
  , _clones(new TClonesArray("TowerInfov4", source.size()))
  , _detector(source.get_detectorid())
{
  for (unsigned int i = 0; i < source.size(); ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerv4::~TowerInfoContainerv4()
{
  delete _clones;
}

void TowerInfoContainerv4::identify(std::ostream& os) const
{
  os << "TowerInfoContainerv4 of size " << size() << std::endl;
}

void TowerInfoContainerv4::Reset()
{
  // clear content of towers in the container for the next event

  for (Int_t i = 0; i < _clones->GetEntriesFast(); ++i)
  {
    TowerInfo *twr = (TowerInfov4 *) _clones->UncheckedAt(i);
    
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

TowerInfov4* TowerInfoContainerv4::get_tower_at_channel(int pos)
{
  return (TowerInfov4*) _clones->At(pos);
}

TowerInfov4* TowerInfoContainerv4::get_tower_at_key(int pos)
{
  int index = decode_key(pos);
  return (TowerInfov4*) _clones->At(index);
}
