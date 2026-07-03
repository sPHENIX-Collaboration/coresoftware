#include "TowerInfoContainerv5.h"
#include "TowerInfov5.h"

#include <TClonesArray.h>
#include <TSystem.h>

#include <cstdlib>

TowerInfoContainerv5::TowerInfoContainerv5(DETECTOR detec)
  : _detector(detec)
{
  int nchannels = get_channels(detec);
  _clones = new TClonesArray("TowerInfov5", nchannels);
  for (int i = 0; i < nchannels; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerv5::TowerInfoContainerv5(const TowerInfoContainerv5& source)
  : TowerInfoContainer(source)
  , _clones(new TClonesArray("TowerInfov5", source.size()))
  , _detector(source.get_detectorid())
{
  for (unsigned int i = 0; i < source.size(); ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

TowerInfoContainerv5::~TowerInfoContainerv5()
{
  delete _clones;
}

void TowerInfoContainerv5::identify(std::ostream& os) const
{
  os << "TowerInfoContainerv5 of size " << size() << std::endl;
}

void TowerInfoContainerv5::Reset()
{
  // clear content of towers in the container for the next event

  for (Int_t i = 0; i < _clones->GetEntriesFast(); ++i)
  {
    TowerInfo *twr = (TowerInfov5 *) _clones->UncheckedAt(i);
    
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

TowerInfov5* TowerInfoContainerv5::get_tower_at_channel(int pos)
{
  return (TowerInfov5*) _clones->At(pos);
}

TowerInfov5* TowerInfoContainerv5::get_tower_at_key(int pos)
{
  int index = decode_key(pos);
  return (TowerInfov5*) _clones->At(index);
}
