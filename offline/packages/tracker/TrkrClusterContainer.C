#include "TrkrClusterContainer.h"
#include "TrkrClusterv1.h"

#include <cstdlib>

using namespace std;

TrkrClusterContainer::TrkrClusterContainer()
{}

void
TrkrClusterContainer::Reset()
{
   while(clusmap.begin() != clusmap.end())
     {
       delete clusmap.begin()->second;
       clusmap.erase(clusmap.begin());
     }
  return;
}

void
TrkrClusterContainer::identify(ostream& os) const
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

TrkrClusterContainer::ConstIterator
TrkrClusterContainer::AddCluster(TrkrCluster *newclus)
{
  TrkrDefs::cluskey key = newclus->get_id();
  if (clusmap.find(key) != clusmap.end())
    {
      TrkrDefUtil util;
      cout << "overwriting clus 0x" << hex << key << dec << endl;
      cout << "tracker ID: " << util.get_trackerid(key) << endl;
    }
  clusmap[key] = newclus;
  return clusmap.find(key);
}

TrkrClusterContainer::ConstIterator
TrkrClusterContainer::AddClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster *newclus)
{
  if(clusmap.find(key)!=clusmap.end())
   {
     cout << "TrkrClusterContainer::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << endl;
     exit(1);
   }
  newclus->set_id(key);
  clusmap[key] = newclus;
  return clusmap.find(key);
}

TrkrClusterContainer::ConstRange 
TrkrClusterContainer::getClusters(const TrkrDefs::TRKRID trackerid) const
{
  // TrkrDefs::cluskey tmp = trackerid;
  // TrkrDefs::cluskey keylow = tmp << TrackerDefs::bitshift_trackerid;
  // TrkrDefs::cluskey keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
//   cout << "keylow: 0x" << hex << keylow << dec << endl;
//   cout << "keyup: 0x" << hex << keyup << dec << endl;

  TrkrDefUtil util;
  TrkrDefs::cluskey keylo = util.get_cluskeylo(trackerid);
  TrkrDefs::cluskey keyhi = util.get_cluskeyhi(trackerid);

  ConstRange retpair;
  retpair.first = clusmap.lower_bound(keylo);
  retpair.second = clusmap.upper_bound(keyhi);
  return retpair;
}

TrkrClusterContainer::ConstRange 
TrkrClusterContainer::getClusters(const TrkrDefs::TRKRID trackerid, const char layer) const
{
  TrkrDefUtil util;
  TrkrDefs::cluskey keylo = util.get_cluskeylo(trackerid, layer);
  TrkrDefs::cluskey keyhi = util.get_cluskeyhi(trackerid, layer);

  ConstRange retpair;
  retpair.first = clusmap.lower_bound(keylo);
  retpair.second = clusmap.upper_bound(keyhi);
  return retpair;
}

TrkrClusterContainer::ConstRange 
TrkrClusterContainer::getClusters( void ) const
{ return std::make_pair( clusmap.begin(), clusmap.end() ); }


TrkrClusterContainer::Iterator 
TrkrClusterContainer::findOrAddCluster(TrkrDefs::cluskey key)
{
  TrkrClusterContainer::Iterator it = clusmap.find(key);
  if(it == clusmap.end())
  {
    clusmap[key] = new TrkrClusterv1();
    it = clusmap.find(key);
    TrkrCluster* mclus = it->second;
    mclus->set_id(key);
  }
  return it;
}

TrkrCluster* 
TrkrClusterContainer::findCluster(TrkrDefs::cluskey key)
{
  TrkrClusterContainer::ConstIterator it = clusmap.find(key);

  if(it != clusmap.end())
    {
      return it->second;
    }

  return 0;
}

