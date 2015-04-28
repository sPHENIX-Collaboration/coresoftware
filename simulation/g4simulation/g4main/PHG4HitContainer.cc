#include "PHG4HitContainer.h"
#include "PHG4Hit.h"
#include "PHG4Hitv1.h"
#include "PHG4HitDefs.h"

#include <phool/phool.h>

#include <cstdlib>

using namespace std;

PHG4HitContainer::PHG4HitContainer()
{
}

void
PHG4HitContainer::Reset()
{
   while(hitmap.begin() != hitmap.end())
     {
       delete hitmap.begin()->second;
       hitmap.erase(hitmap.begin());
     }
  return;
}

void
PHG4HitContainer::identify(ostream& os) const
{
   ConstIterator iter;
   os << "Number of hits: " << size() << endl;
   for (iter = hitmap.begin(); iter != hitmap.end(); iter++)
     {
       os << "hit key 0x" << hex << iter->first << dec << endl;
       (iter->second)->identify();
     }
   set<int>::const_iterator siter;
   os << "Number of layers: " << num_layers() << endl;
   for (siter = layers.begin(); siter != layers.end(); siter++)
     {
       os << "layer : " << *siter << endl;
     }
  return;
}

unsigned long long 
PHG4HitContainer::getmaxkey(const int detid)
{
  ConstRange miter = getHits(detid);
  // first handle 2 special cases where there is no hit in the current layer
  // no hits in this layer and higher layers (lower layers can contain hits)
    if (miter.first == hitmap.end())
    {
      return 0;
    }
  // no hits in this layer - but hits in higher layers
  if (miter.first == miter.second)
    {
      return  0;
    }
  int shiftval = detid << phg4hitdefs::hit_idbits;
  ConstIterator lastlayerentry = miter.second;
  lastlayerentry --;
  unsigned int iret = lastlayerentry->first - shiftval; // subtract layer mask
  return iret;
}


unsigned long long
PHG4HitContainer::genkey(const int detid)
{
  if ((detid >> phg4hitdefs::keybits) > 0)
    {
      cout << " detector id too large: " << detid << endl;
      exit(1);
    }
  unsigned int shiftval = detid << phg4hitdefs::hit_idbits;
  //  cout << "max index: " << (detminmax->second)->first << endl;
  // after removing hits with no energy deposition, we have holes
  // in our hit ranges. This construct will get us the last hit in
  // a layer and return it's hit id. Adding 1 will put us at the end of this layer
  unsigned int hitid = getmaxkey(detid);
  hitid++;
  unsigned int newkey = hitid | shiftval;
  if (hitmap.find(newkey) != hitmap.end())
    {
      cout << PHWHERE << " duplicate key: 0x" 
           << hex << newkey << dec 
	   << " for detector " << detid 
	   << " hitmap.size: " << hitmap.size()
	   << " hitid: " << hitid << " exiting now" << endl;
            exit(1);
    }
  return newkey;
}

std::map<unsigned long long,PHG4Hit *>::const_iterator
PHG4HitContainer::AddHit(const int detid, PHG4Hit *newhit)
{
  unsigned int key = genkey(detid);
  layers.insert(detid);
  newhit->set_hit_id(key);
  hitmap[key] = newhit;
  return hitmap.find(key);
}

PHG4HitContainer::ConstRange PHG4HitContainer::getHits(const int detid) const
{
  if ((detid >> phg4hitdefs::keybits) > 0)
    {
      cout << " detector id too large: " << detid << endl;
      exit(1);
    }
  unsigned int keylow = detid << phg4hitdefs::hit_idbits;
  unsigned int keyup = ((detid + 1)<< phg4hitdefs::hit_idbits) -1 ;
  ConstRange retpair;
  retpair.first = hitmap.lower_bound(keylow);
  retpair.second = hitmap.upper_bound(keyup);
  return retpair;
}

PHG4HitContainer::ConstRange PHG4HitContainer::getHits( void ) const
{ return std::make_pair( hitmap.begin(), hitmap.end() ); }


PHG4HitContainer::Iterator PHG4HitContainer::findOrAddHit(unsigned long long key)
{
  PHG4HitContainer::Iterator it = hitmap.find(key);
  if(it == hitmap.end())
  {
    hitmap[key] = new PHG4Hitv1();
    it = hitmap.find(key);
    PHG4Hit* mhit = it->second;
    mhit->set_hit_id(key);
    mhit->set_edep(0.);
    layers.insert(mhit->get_layer()); // add layer to our set of layers
  }
  return it;
}

PHG4Hit* PHG4HitContainer::findHit(unsigned long long key)
{
  PHG4HitContainer::ConstIterator it = hitmap.find(key);
  if(it != hitmap.end())
    {
      return it->second;
    }
    
  return 0;
}

void
PHG4HitContainer::RemoveZeroEDep()
{
  //  unsigned int hitsbef = hitmap.size();
  Iterator itr = hitmap.begin();
  Iterator last = hitmap.end();
  for (; itr != last; )
    {
      PHG4Hit *hit = itr->second;
      if (hit->get_edep() == 0)
        {
          delete hit;
          hitmap.erase(itr++);
        }
      else
        {
          ++itr;
        }
    }
//   unsigned int hitsafter = hitmap.size();
//   cout << "hist before: " << hitsbef
//        << ", hits after: " << hitsafter << endl;
  return;
}

