#include "TrackerHitContainer.h"
#include "TrackerHitv1.h"

#include <cstdlib>

using namespace std;

TrackerHitContainer::TrackerHitContainer()
{}

void
TrackerHitContainer::Reset()
{
   while(hitmap.begin() != hitmap.end())
     {
       delete hitmap.begin()->second;
       hitmap.erase(hitmap.begin());
     }
  return;
}

void
TrackerHitContainer::identify(ostream& os) const
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

TrackerHitContainer::ConstIterator
TrackerHitContainer::AddHit(TrackerHit *newhit)
{
  TrackerDefs::keytype key = newhit->get_hitid();
  if (hitmap.find(key) != hitmap.end())
    {
      cout << "overwriting hit 0x" << hex << key << dec << endl;
      cout << "tracker id: " << TrackerDefs::get_trackerid(key) << endl;
    }
  hitmap[key] = newhit;
  return hitmap.find(key);
}

TrackerHitContainer::ConstIterator
TrackerHitContainer::AddHitSpecifyKey(const TrackerDefs::keytype key, TrackerHit *newhit)
{
  if(hitmap.find(key)!=hitmap.end())
   {
     cout << "TrackerHitContainer::AddHitSpecifyKey: duplicate key: " << key << " exiting now" << endl;
     exit(1);
   }
  newhit->set_hitid(key);
  hitmap[key] = newhit;
  return hitmap.find(key);
}

TrackerHitContainer::ConstRange 
TrackerHitContainer::getHits(const TrackerDefs::TRACKERID trackerid) const
{
  TrackerDefs::keytype tmp = trackerid;
  TrackerDefs::keytype keylow = tmp << TrackerDefs::bitshift_trackerid;
  TrackerDefs::keytype keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
//   cout << "keylow: 0x" << hex << keylow << dec << endl;
//   cout << "keyup: 0x" << hex << keyup << dec << endl;
  ConstRange retpair;
  retpair.first = hitmap.lower_bound(keylow);
  retpair.second = hitmap.upper_bound(keyup);
  return retpair;
}

TrackerHitContainer::ConstRange 
TrackerHitContainer::getHits(const TrackerDefs::TRACKERID trackerid, 
  const char layer) const
{
  TrackerDefs::keytype tmp = trackerid;
  TrackerDefs::keytype keylow = (tmp << TrackerDefs::bitshift_trackerid);
  tmp = layer;
  keylow |= (tmp << TrackerDefs::bitshift_layer);
  TrackerDefs::keytype keyup = ((tmp + 1)<< TrackerDefs::bitshift_layer) -1 ;
//   cout << "keylow: 0x" << hex << keylow << dec << endl;
//   cout << "keyup: 0x" << hex << keyup << dec << endl;
  ConstRange retpair;
  retpair.first = hitmap.lower_bound(keylow);
  retpair.second = hitmap.upper_bound(keyup);
  return retpair;
}

TrackerHitContainer::ConstRange 
TrackerHitContainer::getHits( void ) const
{ return std::make_pair( hitmap.begin(), hitmap.end() ); }


TrackerHitContainer::Iterator 
TrackerHitContainer::findOrAddHit(TrackerDefs::keytype key)
{
  TrackerHitContainer::Iterator it = hitmap.find(key);
  if(it == hitmap.end())
  {
    hitmap[key] = new TrackerHitv1();
    it = hitmap.find(key);
    TrackerHit* mhit = it->second;
    mhit->set_hitid(key);
  }
  return it;
}

TrackerHit* 
TrackerHitContainer::findHit(TrackerDefs::keytype key)
{
  TrackerHitContainer::ConstIterator it = hitmap.find(key);

  if(it != hitmap.end())
    {
      return it->second;
    }

  return 0;
}
