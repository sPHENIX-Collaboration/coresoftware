#include "PHG4HitContainer.h"

#include "PHG4Hit.h"
#include "PHG4Hitv1.h"

#include <phool/phool.h>

#include <TSystem.h>

#include <cstdlib>

using namespace std;

PHG4HitContainer::PHG4HitContainer()
  : id(-1), hitmap(), layers()
{
}

PHG4HitContainer::PHG4HitContainer(const std::string &nodename)
  : id(PHG4HitDefs::get_volume_id(nodename)), hitmap(), layers()
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
   for (iter = hitmap.begin(); iter != hitmap.end(); ++iter)
     {
       os << "hit key 0x" << hex << iter->first << dec << endl;
       (iter->second)->identify();
     }
   set<unsigned int>::const_iterator siter;
   os << "Number of layers: " << num_layers() << endl;
   for (siter = layers.begin(); siter != layers.end(); ++siter)
     {
       os << "layer : " << *siter << endl;
     }
  return;
}

PHG4HitDefs::keytype
PHG4HitContainer::getmaxkey(const unsigned int detid)
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
  PHG4HitDefs::keytype detidlong = detid;
  PHG4HitDefs::keytype shiftval = detidlong << PHG4HitDefs::hit_idbits;
  ConstIterator lastlayerentry = miter.second;
  --lastlayerentry;
  PHG4HitDefs::keytype iret = lastlayerentry->first - shiftval; // subtract layer mask
  return iret;
}


PHG4HitDefs::keytype
PHG4HitContainer::genkey(const unsigned int detid)
{
  PHG4HitDefs::keytype detidlong = detid;
  if ((detidlong >> PHG4HitDefs::keybits) > 0)
    {
      cout << PHWHERE << " detector id too large: " << detid << endl;
      gSystem->Exit(1);
    }
  PHG4HitDefs::keytype shiftval = detidlong << PHG4HitDefs::hit_idbits;
  //  cout << "max index: " << (detminmax->second)->first << endl;
  // after removing hits with no energy deposition, we have holes
  // in our hit ranges. This construct will get us the last hit in
  // a layer and return it's hit id. Adding 1 will put us at the end of this layer
  PHG4HitDefs::keytype hitid = getmaxkey(detid);
  hitid++;
  PHG4HitDefs::keytype newkey = hitid | shiftval;
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

PHG4HitContainer::ConstIterator
PHG4HitContainer::AddHit(PHG4Hit *newhit)
{
  PHG4HitDefs::keytype key = newhit->get_hit_id();
  if (hitmap.find(key) != hitmap.end())
    {
      cout << "hit with id  0x" << hex << key << dec << " exists already" << endl;
      return hitmap.find(key);
    }
  PHG4HitDefs::keytype detidlong = key >>  PHG4HitDefs::hit_idbits;
  unsigned int detid = detidlong;
  layers.insert(detid);
  return hitmap.insert( std::make_pair( key, newhit ) ).first;
}

PHG4HitContainer::ConstIterator
PHG4HitContainer::AddHit(const unsigned int detid, PHG4Hit *newhit)
{
  PHG4HitDefs::keytype key = genkey(detid);
  layers.insert(detid);
  newhit->set_hit_id(key);
  return hitmap.insert( std::make_pair( key, newhit ) ).first;
}

PHG4HitContainer::ConstRange PHG4HitContainer::getHits(const unsigned int detid) const
{
  PHG4HitDefs::keytype detidlong = detid;
  if ((detidlong >> PHG4HitDefs::keybits) > 0)
    {
      cout << " detector id too large: " << detid << endl;
      exit(1);
    }
  PHG4HitDefs::keytype keylow = detidlong << PHG4HitDefs::hit_idbits;
  PHG4HitDefs::keytype keyup = ((detidlong + 1) << PHG4HitDefs::hit_idbits) -1 ;
  ConstRange retpair;
  retpair.first = hitmap.lower_bound(keylow);
  retpair.second = hitmap.upper_bound(keyup);
  return retpair;
}

PHG4HitContainer::ConstRange PHG4HitContainer::getHits( void ) const
{ return std::make_pair( hitmap.begin(), hitmap.end() ); }


PHG4HitContainer::Iterator PHG4HitContainer::findOrAddHit(PHG4HitDefs::keytype key)
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

PHG4Hit* PHG4HitContainer::findHit(PHG4HitDefs::keytype key)
{
  PHG4HitContainer::ConstIterator it = hitmap.find(key);
  if(it != hitmap.end())
    {
      return it->second;
    }
    
  return nullptr;
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

