#include "TrkrHitSetContainer.h"
#include "TrkrHitSetv1.h"
#include "TrkrDefUtil.h"

#include <cstdlib>

using namespace std;

TrkrHitSetContainer::TrkrHitSetContainer()
{}

void
TrkrHitSetContainer::Reset()
{
   while(hitmap.begin() != hitmap.end())
     {
       delete hitmap.begin()->second;
       hitmap.erase(hitmap.begin());
     }
  return;
}

void
TrkrHitSetContainer::identify(ostream& os) const
{
   ConstIterator iter;
   os << "Number of hits: " << size() << endl;
   for (iter = hitmap.begin(); iter != hitmap.end(); ++iter)
     {
       os << "hit key 0x" << hex << iter->first << dec << endl;
       (iter->second)->identify();
     }
  return;
}

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::AddHitSet(TrkrHitSet *newhit)
{
  TrkrDefs::hitsetkey key = newhit->get_hitid();
  if (hitmap.find(key) != hitmap.end())
    {
      TrkrDefUtil util;
      cout << "overwriting hit 0x" << hex << key << dec << endl;
      cout << "tracker id: " << util.get_trackerid(key) << endl;
    }
  hitmap[key] = newhit;
  return hitmap.find(key);
}

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::AddHitSetSpecifyKey(const TrkrDefs::hitsetkey key, TrkrHitSet *newhit)
{
  if(hitmap.find(key)!=hitmap.end())
   {
     cout << "TrkrHitSetContainer::AddHitSpecifyKey: duplicate key: " << key << " exiting now" << endl;
     exit(1);
   }
  newhit->set_hitid(key);
  hitmap[key] = newhit;
  return hitmap.find(key);
}

TrkrHitSetContainer::ConstRange 
TrkrHitSetContainer::getHitSets(const TrkrDefs::TRKRID trackerid) const
{

  // TrkrDefs::hitsetkey tmp = trackerid;
  // TrkrDefs::hitsetkey keylow = tmp << TrackerDefs::bitshift_trackerid;
  // TrkrDefs::hitsetkey keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
//   cout << "keylow: 0x" << hex << keylow << dec << endl;
//   cout << "keyup: 0x" << hex << keyup << dec << endl;

  TrkrDefUtil util;
  TrkrDefs::hitsetkey keylo = util.get_hitsetkeylo(trackerid);
  TrkrDefs::hitsetkey keyhi = util.get_hitsetkeyhi(trackerid);

  ConstRange retpair;
  retpair.first = hitmap.lower_bound(keylo);
  retpair.second = hitmap.upper_bound(keyhi);
  return retpair;
}

TrkrHitSetContainer::ConstRange 
TrkrHitSetContainer::getHitSets(const TrkrDefs::TRKRID trackerid, 
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
  TrkrDefs::hitsetkey keylo = util.get_hitsetkeylo(trackerid, layer);
  TrkrDefs::hitsetkey keyhi = util.get_hitsetkeyhi(trackerid, layer);

  ConstRange retpair;
  retpair.first = hitmap.lower_bound(keylo);
  retpair.second = hitmap.upper_bound(keyhi);
  return retpair;
}

TrkrHitSetContainer::ConstRange 
TrkrHitSetContainer::getHitSets( void ) const
{ return std::make_pair( hitmap.begin(), hitmap.end() ); }


TrkrHitSetContainer::Iterator 
TrkrHitSetContainer::findOrAddHitSet(TrkrDefs::hitsetkey key)
{
  TrkrHitSetContainer::Iterator it = hitmap.find(key);
  if(it == hitmap.end())
  {
    hitmap[key] = new TrkrHitSetv1();
    it = hitmap.find(key);
    TrkrHitSet* mhit = it->second;
    mhit->set_hitid(key);
  }
  return it;
}

TrkrHitSet* 
TrkrHitSetContainer::findHitSet(TrkrDefs::hitsetkey key)
{
  TrkrHitSetContainer::ConstIterator it = hitmap.find(key);

  if(it != hitmap.end())
    {
      return it->second;
    }

  return 0;
}
