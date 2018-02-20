#include "TrkrClusterContainer.h"
#include "TrkrClusterv1.h"

#include <cstdlib>

using namespace std;

TrkrClusterContainer::TrkrClusterContainer()
{
}

void TrkrClusterContainer::Reset()
{
  while (clusmap_.begin() != clusmap_.end())
  {
    delete clusmap_.begin()->second;
    clusmap_.erase(clusmap_.begin());
  }
  return;
}

void TrkrClusterContainer::identify(ostream& os) const
{
  ConstIterator iter;
  os << "Number of clusters: " << size() << endl;
  for (iter = clusmap_.begin(); iter != clusmap_.end(); ++iter)
  {
    os << "clus key 0x" << hex << iter->first << dec << endl;
    (iter->second)->identify();
  }
  return;
}

TrkrClusterContainer::ConstIterator
TrkrClusterContainer::AddCluster(TrkrCluster* newclus)
{
  TrkrDefs::cluskey key = newclus->GetClusKey();
  if (clusmap_.find(key) != clusmap_.end())
  {
    TrkrDefUtil util;
    cout << "overwriting clus 0x" << hex << key << dec << endl;
    cout << "tracker ID: " << util.GetTrkrId(key) << endl;
  }
  clusmap_[key] = newclus;
  return clusmap_.find(key);
}

TrkrClusterContainer::ConstIterator
TrkrClusterContainer::AddClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster* newclus)
{
  if (clusmap_.find(key) != clusmap_.end())
  {
    cout << "TrkrClusterContainer::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << endl;
    exit(1);
  }
  newclus->SetClusKey(key);
  clusmap_[key] = newclus;
  return clusmap_.find(key);
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::GetClusters(const TrkrDefs::TRKRID trackerid) const
{
  // TrkrDefs::cluskey tmp = trackerid;
  // TrkrDefs::cluskey keylow = tmp << TrackerDefs::bitshift_trackerid;
  // TrkrDefs::cluskey keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
  //   cout << "keylow: 0x" << hex << keylow << dec << endl;
  //   cout << "keyup: 0x" << hex << keyup << dec << endl;

  TrkrDefUtil util;
  TrkrDefs::cluskey keylo = util.GetClusKeyLo(trackerid);
  TrkrDefs::cluskey keyhi = util.GetClusKeyHi(trackerid);

  ConstRange retpair;
  retpair.first = clusmap_.lower_bound(keylo);
  retpair.second = clusmap_.upper_bound(keyhi);
  return retpair;
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::GetClusters(const TrkrDefs::TRKRID trackerid, const char layer) const
{
  TrkrDefUtil util;
  TrkrDefs::cluskey keylo = util.GetClusKeyLo(trackerid, layer);
  TrkrDefs::cluskey keyhi = util.GetClusKeyHi(trackerid, layer);

  ConstRange retpair;
  retpair.first = clusmap_.lower_bound(keylo);
  retpair.second = clusmap_.upper_bound(keyhi);
  return retpair;
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::GetClusters(void) const
{
  return std::make_pair(clusmap_.begin(), clusmap_.end());
}

TrkrClusterContainer::Iterator
TrkrClusterContainer::FindOrAddCluster(TrkrDefs::cluskey key)
{
  TrkrClusterContainer::Iterator it = clusmap_.find(key);
  if (it == clusmap_.end())
  {
    clusmap_[key] = new TrkrClusterv1();
    it = clusmap_.find(key);
    TrkrCluster* mclus = it->second;
    mclus->SetClusKey(key);
  }
  return it;
}

TrkrCluster*
TrkrClusterContainer::FindCluster(TrkrDefs::cluskey key)
{
  TrkrClusterContainer::ConstIterator it = clusmap_.find(key);

  if (it != clusmap_.end())
  {
    return it->second;
  }

  return 0;
}
