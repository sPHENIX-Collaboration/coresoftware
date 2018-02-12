#include "TrkrHitSetContainer.h"
#include "TrkrHitSetv1.h"

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
TrkrHitSetContainer::AddHit(TrkrHitSet *newhit)
{
  TrackerDefs::hitkeytype key = newhit->get_hitid();
  if (hitmap.find(key) != hitmap.end())
    {
      cout << "overwriting hit 0x" << hex << key << dec << endl;
      cout << "tracker id: " << TrackerDefs::get_trackerid(key) << endl;
    }
  hitmap[key] = newhit;
  return hitmap.find(key);
}

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::AddHitSpecifyKey(const TrackerDefs::hitkeytype key, TrkrHitSet *newhit)
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
TrkrHitSetContainer::getHits(const TrackerDefs::TRACKERID trackerid) const
{
  TrackerDefs::hitkeytype tmp = trackerid;
  TrackerDefs::hitkeytype keylow = tmp << TrackerDefs::bitshift_trackerid;
  TrackerDefs::hitkeytype keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
//   cout << "keylow: 0x" << hex << keylow << dec << endl;
//   cout << "keyup: 0x" << hex << keyup << dec << endl;
  ConstRange retpair;
  retpair.first = hitmap.lower_bound(keylow);
  retpair.second = hitmap.upper_bound(keyup);
  return retpair;
}

TrkrHitSetContainer::ConstRange 
TrkrHitSetContainer::getHits(const TrackerDefs::TRACKERID trackerid, 
  const char layer) const
{
  TrackerDefs::hitkeytype tmp = trackerid;
  TrackerDefs::hitkeytype keylow = (tmp << TrackerDefs::bitshift_trackerid);
  tmp = layer;
  keylow |= (tmp << TrackerDefs::bitshift_layer);
  TrackerDefs::hitkeytype keyup = ((tmp + 1)<< TrackerDefs::bitshift_layer) -1 ;
//   cout << "keylow: 0x" << hex << keylow << dec << endl;
//   cout << "keyup: 0x" << hex << keyup << dec << endl;
  ConstRange retpair;
  retpair.first = hitmap.lower_bound(keylow);
  retpair.second = hitmap.upper_bound(keyup);
  return retpair;
}

TrkrHitSetContainer::ConstRange 
TrkrHitSetContainer::getHits( void ) const
{ return std::make_pair( hitmap.begin(), hitmap.end() ); }


TrkrHitSetContainer::Iterator 
TrkrHitSetContainer::findOrAddHit(TrackerDefs::hitkeytype key)
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
TrkrHitSetContainer::findHit(TrackerDefs::hitkeytype key)
{
  TrkrHitSetContainer::ConstIterator it = hitmap.find(key);

  if(it != hitmap.end())
    {
      return it->second;
    }

  return 0;
}
