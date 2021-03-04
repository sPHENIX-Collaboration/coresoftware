/**
 * @file trackbase/TrkrClusterContainer.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrClusterContainer
 */
#include "TrkrClusterContainer.h"
#include "TrkrCluster.h"
#include "TrkrClusterv1.h"
#include "TrkrDefs.h"

#include <tpc/TpcDefs.h>

//#include <boost/range.hpp>
//#include <boost/range/join.hpp>
#include <cstdlib>

void TrkrClusterContainer::Reset()
{
  for(unsigned int layer = 0;layer < max_layer; layer++){
    for(unsigned int phi_segment = 0;phi_segment < max_phisegment;phi_segment++){
      for(unsigned int z_segment = 0; z_segment < max_zsegment; z_segment++){
	while (m_clusmap[layer][phi_segment][z_segment].begin() != m_clusmap[layer][phi_segment][z_segment].end())
	  {
	    delete m_clusmap[layer][phi_segment][z_segment].begin()->second;
	    m_clusmap[layer][phi_segment][z_segment].erase(m_clusmap[layer][phi_segment][z_segment].begin());
	  }
      }
    }
  }
  return;
}

void TrkrClusterContainer::identify(std::ostream& os) const
{
  os << "-----TrkrClusterContainer-----" << std::endl;
  ConstIterator iter;
  os << "Number of clusters: " << size() << std::endl;
  for(unsigned int layer = 0;layer < max_layer; layer++){
    for(unsigned int phi_segment = 0;phi_segment < max_phisegment;phi_segment++){
      for(unsigned int z_segment = 0; z_segment < max_zsegment; z_segment++){
	for (iter = m_clusmap[layer][phi_segment][z_segment].begin(); 
	     iter != m_clusmap[layer][phi_segment][z_segment].end();
	     ++iter)
	  {
	    int layer = TrkrDefs::getLayer(iter->first);
	    os << "clus key " << iter->first  << " layer " << layer << std::endl;
	    (iter->second)->identify();
	  }
      }
    }
  }
  os << "------------------------------" << std::endl;
  return;
}

void TrkrClusterContainer::removeCluster(TrkrDefs::cluskey key){
  unsigned int layer = TrkrDefs::getLayer(key);
  unsigned int sector = TpcDefs::getSectorId(key);
  unsigned int side = TpcDefs::getSide(key); 
  m_clusmap[layer][sector][side].erase(key);
}

void TrkrClusterContainer::removeCluster(TrkrCluster *clus){
  TrkrDefs::cluskey key = clus->getClusKey();
  unsigned int layer = TrkrDefs::getLayer(key);
  unsigned int sector = TpcDefs::getSectorId(key);
  unsigned int side = TpcDefs::getSide(key); 
  m_clusmap[layer][sector][side].erase(key);

}
 
TrkrClusterContainer::ConstIterator
TrkrClusterContainer::addCluster(TrkrCluster* newclus)
{
  return addClusterSpecifyKey(newclus->getClusKey(), newclus);
}

TrkrClusterContainer::ConstIterator
TrkrClusterContainer::addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster* newclus)
{
  unsigned int layer = TrkrDefs::getLayer(key);
  unsigned int sector = TpcDefs::getSectorId(key);
  unsigned int side = TpcDefs::getSide(key); 
  auto ret = m_clusmap[layer][sector][side].insert(std::make_pair(key, newclus));
  if ( !ret.second )
  {
    std::cout << "TrkrClusterContainer::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
  else
  {
    return ret.first;
  }
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::getClusters(const TrkrDefs::TrkrId trackerid) const
{
  // find first layer of this trackerId
  unsigned int min_layer = UINT_MAX;

  for(unsigned int layer = 0; layer < max_layer; layer++){
    for(unsigned int phi_segment = 0;phi_segment < max_phisegment;phi_segment++){
      for(unsigned int z_segment = 0; z_segment < max_zsegment; z_segment++){
	if(m_clusmap[layer][phi_segment][z_segment].begin()!=
	   m_clusmap[layer][phi_segment][z_segment].end()){
	  auto iter = m_clusmap[layer][phi_segment][z_segment].begin();
	  TrkrDefs::TrkrId this_trackerid = (TrkrDefs::TrkrId) TrkrDefs::getTrkrId(iter->first);
	  if(this_trackerid == trackerid)
	    min_layer = TrkrDefs::getLayer(iter->first);
	}
      }
    }
  }
  // TrkrDefs::cluskey tmp = trackerid;
  // TrkrDefs::cluskey keylow = tmp << TrackerDefs::bitshift_trackerid;
  // TrkrDefs::cluskey keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
  //   std::cout << "keylow: 0x" << hex << keylow << dec << std::endl;
  //   std::cout << "keyup: 0x" << hex << keyup << dec << std::endl;

  TrkrDefs::cluskey keylo = TrkrDefs::getClusKeyLo(trackerid);
  TrkrDefs::cluskey keyhi = TrkrDefs::getClusKeyHi(trackerid);

  ConstRange retpair;
  if(min_layer != UINT_MAX){
    retpair.first  = m_clusmap[min_layer][0][0].lower_bound(keylo);
    retpair.second = m_clusmap[min_layer][0][0].upper_bound(keyhi);
  }
  return retpair;
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::getClusters(const TrkrDefs::TrkrId trackerid, const unsigned int layer) const
{
  TrkrDefs::cluskey keylo = TrkrDefs::getClusKeyLo(trackerid, layer);
  TrkrDefs::cluskey keyhi = TrkrDefs::getClusKeyHi(trackerid, layer);

  ConstRange retpair;

  retpair.first = m_clusmap[layer][0][0].lower_bound(keylo);
  retpair.second = m_clusmap[layer][0][0].upper_bound(keyhi);

  return retpair;
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::getClusters(const unsigned int layer) const
{
  // TrkrDefs::cluskey keylo = TrkrDefs::getClusKeyLo(trackerid, layer);
  // TrkrDefs::cluskey keyhi = TrkrDefs::getClusKeyHi(trackerid, layer);
  // retpair.first = m_clusmap.lower_bound(keylo);
  // retpair.second = m_clusmap.upper_bound(keyhi);

  ConstRange retpair;
  retpair.first  = m_clusmap[layer][0][0].begin();
  retpair.second = m_clusmap[layer][0][0].end();


  return retpair;
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::getClusters(void) const
{
  unsigned int layer = 7;
  return std::make_pair(m_clusmap[layer][0][0].begin(), m_clusmap[layer][0][0].end());
}

TrkrClusterContainer::Iterator
TrkrClusterContainer::findOrAddCluster(TrkrDefs::cluskey key){
  unsigned int layer  = TrkrDefs::getLayer(key);
  unsigned int sector = TpcDefs::getSectorId(key);
  unsigned int side   = TpcDefs::getSide(key); 
  TrkrClusterContainer::Iterator it = m_clusmap[layer][sector][side].find(key);
  if (it == m_clusmap[layer][sector][side].end())
  {
    // add new cluster and set its key
    auto ret = m_clusmap[layer][sector][side].insert(std::make_pair(key, new TrkrClusterv1()));
    (ret.first->second)->setClusKey(key);
    it = ret.first;
  }
  return it;
}

TrkrCluster* TrkrClusterContainer::findCluster(TrkrDefs::cluskey key){
  unsigned int layer  = TrkrDefs::getLayer(key);
  unsigned int sector = TpcDefs::getSectorId(key);
  unsigned int side   = TpcDefs::getSide(key); 
  TrkrClusterContainer::Iterator it = m_clusmap[layer][sector][side].find(key);

  if (it != m_clusmap[layer][sector][side].end())
  {
    return it->second;
  }

  return 0;
}
