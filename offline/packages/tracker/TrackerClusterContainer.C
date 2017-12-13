#include "TrackerClusterContainer.h"
#include "TrackerClusterv1.h"
#include "TrackerDefs.h"

#include <cstdlib>

using namespace std;

TrackerClusterContainer::TrackerClusterContainer()
{}

void
TrackerClusterContainer::Reset()
{
   while(clusmap.begin() != clusmap.end())
     {
       delete clusmap.begin()->second;
       clusmap.erase(clusmap.begin());
     }
  return;
}

void
TrackerClusterContainer::identify(ostream& os) const
{
   ConstIterator iter;
   os << "Number of clusters: " << size() << endl;
   for (iter = clusmap.begin(); iter != clusmap.end(); ++iter)
     {
       os << "clus key 0x" << hex << iter->first << dec << endl;
       (iter->second)->identify();
     }
  return;
}

TrackerClusterContainer::ConstIterator
TrackerClusterContainer::AddCluster(TrackerCluster *newclus)
{
  TrackerDefs::keytype key = newclus->get_id();
  if (clusmap.find(key) != clusmap.end())
    {
      cout << "overwriting clus 0x" << hex << key << dec << endl;
      cout << "tracker ID: " << TrackerDefs::get_trackerid(key) << endl;
    }
  clusmap[key] = newclus;
  return clusmap.find(key);
}

TrackerClusterContainer::ConstIterator
TrackerClusterContainer::AddClusterSpecifyKey(const TrackerDefs::keytype key, TrackerCluster *newclus)
{
  if(clusmap.find(key)!=clusmap.end())
   {
     cout << "TrackerClusterContainer::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << endl;
     exit(1);
   }
  newclus->set_id(key);
  clusmap[key] = newclus;
  return clusmap.find(key);
}

TrackerClusterContainer::ConstRange 
TrackerClusterContainer::getClusters(const TrackerDefs::TRACKERID trackerid) const
{
  TrackerDefs::keytype tmp = trackerid;
  TrackerDefs::keytype keylow = tmp << TrackerDefs::bitshift_trackerid;
  TrackerDefs::keytype keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
//   cout << "keylow: 0x" << hex << keylow << dec << endl;
//   cout << "keyup: 0x" << hex << keyup << dec << endl;
  ConstRange retpair;
  retpair.first = clusmap.lower_bound(keylow);
  retpair.second = clusmap.upper_bound(keyup);
  return retpair;
}

TrackerClusterContainer::ConstRange 
TrackerClusterContainer::getClusters( void ) const
{ return std::make_pair( clusmap.begin(), clusmap.end() ); }


TrackerClusterContainer::Iterator 
TrackerClusterContainer::findOrAddCluster(TrackerDefs::keytype key)
{
  TrackerClusterContainer::Iterator it = clusmap.find(key);
  if(it == clusmap.end())
  {
    clusmap[key] = new TrackerClusterv1();
    it = clusmap.find(key);
    TrackerCluster* mclus = it->second;
    mclus->set_id(key);
  }
  return it;
}

TrackerCluster* 
TrackerClusterContainer::findCluster(TrackerDefs::keytype key)
{
  TrackerClusterContainer::ConstIterator it = clusmap.find(key);

  if(it != clusmap.end())
    {
      return it->second;
    }

  return 0;
}

