#include "TrkrDefUtil.h"
#include "TrkrHitSetContainer.h"
#include "TrkrHitSetv1.h"

#include <cstdlib>

using namespace std;

TrkrHitSetContainer::TrkrHitSetContainer()
{
}

void TrkrHitSetContainer::Reset()
{
  while (hitmap_.begin() != hitmap_.end())
  {
    delete hitmap_.begin()->second;
    hitmap_.erase(hitmap_.begin());
  }
  return;
}

void TrkrHitSetContainer::identify(ostream& os) const
{
  ConstIterator iter;
  os << "Number of hits: " << size() << endl;
  for (iter = hitmap_.begin(); iter != hitmap_.end(); ++iter)
  {
    os << "hit key 0x" << hex << iter->first << dec << endl;
    (iter->second)->identify();
  }
  return;
}

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::AddHitSet(TrkrHitSet* newhit)
{
  TrkrDefs::hitsetkey key = newhit->GetHitSetKey();
  if (hitmap_.find(key) != hitmap_.end())
  {
    TrkrDefUtil util;
    cout << "overwriting hit 0x" << hex << key << dec << endl;
    cout << "tracker id: " << util.GetTrkrId(key) << endl;
  }
  hitmap_[key] = newhit;
  return hitmap_.find(key);
}

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::AddHitSetSpecifyKey(const TrkrDefs::hitsetkey key, TrkrHitSet* newhit)
{
  if (hitmap_.find(key) != hitmap_.end())
  {
    cout << "TrkrHitSetContainer::AddHitSpecifyKey: duplicate key: " << key << " exiting now" << endl;
    exit(1);
  }
  newhit->SetHitSetKey(key);
  hitmap_[key] = newhit;
  return hitmap_.find(key);
}

TrkrHitSetContainer::ConstRange
TrkrHitSetContainer::GetHitSets(const TrkrDefs::TRKRID trackerid) const
{
  // TrkrDefs::hitsetkey tmp = trackerid;
  // TrkrDefs::hitsetkey keylow = tmp << TrackerDefs::bitshift_trackerid;
  // TrkrDefs::hitsetkey keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
  //   cout << "keylow: 0x" << hex << keylow << dec << endl;
  //   cout << "keyup: 0x" << hex << keyup << dec << endl;

  TrkrDefUtil util;
  TrkrDefs::hitsetkey keylo = util.GetHitSetKeyLo(trackerid);
  TrkrDefs::hitsetkey keyhi = util.GetHitSetKeyHi(trackerid);

  ConstRange retpair;
  retpair.first = hitmap_.lower_bound(keylo);
  retpair.second = hitmap_.upper_bound(keyhi);
  return retpair;
}

TrkrHitSetContainer::ConstRange
TrkrHitSetContainer::GetHitSets(const TrkrDefs::TRKRID trackerid,
                                const char layer) const
{
  // TrkrDefs::hitsetkey tmp = trackerid;
  // TrkrDefs::hitsetkey keylow = (tmp << TrackerDefs::bitshift_trackerid);
  // tmp = layer;
  // keylow |= (tmp << TrackerDefs::bitshift_layer);
  // TrkrDefs::hitsetkey keyup = ((tmp + 1)<< TrackerDefs::bitshift_layer) -1 ;
  //   cout << "keylow: 0x" << hex << keylow << dec << endl;
  //   cout << "keyup: 0x" << hex << keyup << dec << endl;

  TrkrDefUtil util;
  TrkrDefs::hitsetkey keylo = util.GetHitSetKeyLo(trackerid, layer);
  TrkrDefs::hitsetkey keyhi = util.GetHitSetKeyHi(trackerid, layer);

  ConstRange retpair;
  retpair.first = hitmap_.lower_bound(keylo);
  retpair.second = hitmap_.upper_bound(keyhi);
  return retpair;
}

TrkrHitSetContainer::ConstRange
TrkrHitSetContainer::GetHitSets(void) const
{
  return std::make_pair(hitmap_.begin(), hitmap_.end());
}

TrkrHitSetContainer::Iterator
TrkrHitSetContainer::FindOrAddHitSet(TrkrDefs::hitsetkey key)
{
  TrkrHitSetContainer::Iterator it = hitmap_.find(key);
  if (it == hitmap_.end())
  {
    hitmap_[key] = new TrkrHitSetv1();
    it = hitmap_.find(key);
    TrkrHitSet* mhit = it->second;
    mhit->SetHitSetKey(key);
  }
  return it;
}

TrkrHitSet*
TrkrHitSetContainer::FindHitSet(TrkrDefs::hitsetkey key)
{
  TrkrHitSetContainer::ConstIterator it = hitmap_.find(key);

  if (it != hitmap_.end())
  {
    return it->second;
  }

  return 0;
}
