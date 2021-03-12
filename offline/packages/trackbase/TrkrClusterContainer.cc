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

void TrkrClusterContainer::print() const
{
  for(unsigned int layer = 0;layer < max_layer; layer++){
    for(unsigned int phi_segment = 0;phi_segment < max_phisegment;phi_segment++){
      for(unsigned int z_segment = 0; z_segment < max_zsegment; z_segment++){
	auto iter = m_clusmap[layer][phi_segment][z_segment].begin(); 
	int flayer = TrkrDefs::getLayer(iter->first);
	unsigned int fsector = TrkrDefs::getPhiElement(iter->first);
	unsigned int fside = TrkrDefs::getZElement(iter->first); 
	std::cout << "layer: " << layer << " | " << flayer
		  << " phi_seg: " << phi_segment << " | " << fsector
		  << " z_seg: " << z_segment << " | " << fside
		  << " nclu: " << m_clusmap[layer][phi_segment][z_segment].size()
		  << std::endl;

      }
    }
  }
  return;
}

void TrkrClusterContainer::removeCluster(TrkrDefs::cluskey key){
  unsigned int layer = TrkrDefs::getLayer(key);
  unsigned int sector = TrkrDefs::getPhiElement(key);
  unsigned int side = TrkrDefs::getZElement(key); 
  m_clusmap[layer][sector][side].erase(key);
}

void TrkrClusterContainer::removeCluster(TrkrCluster *clus){
  TrkrDefs::cluskey key = clus->getClusKey();
  unsigned int layer = TrkrDefs::getLayer(key);
  unsigned int sector = TrkrDefs::getPhiElement(key);
  unsigned int side = TrkrDefs::getZElement(key); 
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
  unsigned int sector = TrkrDefs::getPhiElement(key);
  unsigned int side = TrkrDefs::getZElement(key); 
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
TrkrClusterContainer::getClusters() const
{
  std::cout << "DEPRECATED METHOD TO ACCESS CLUSTER CONTAINER!!!" << std::endl;
  ConstRange retpair;
  retpair.first  = m_clusmap[7][0][0].begin();
  retpair.second = m_clusmap[7][0][0].end();

  return retpair;
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::getClusters(TrkrDefs::hitsetkey hitsetkey) const
{
  const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  const unsigned int side  = TrkrDefs::getZElement(hitsetkey);
  const unsigned int sector= TrkrDefs::getPhiElement(hitsetkey);
  ConstRange retpair;
  if(sector<max_phisegment&&side<max_zsegment){
    retpair.first  = m_clusmap[layer][sector][side].begin();
    retpair.second = m_clusmap[layer][sector][side].end();
  }else{
    retpair.first  = m_clusmap[layer][sector][side].begin();
    retpair.second = m_clusmap[layer][sector][side].begin();
  }
  return retpair;
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::getClusters(const unsigned int layer, const unsigned int phi_segment, const unsigned int z_segment) const
{
  // TrkrDefs::cluskey keylo = TrkrDefs::getClusKeyLo(trackerid, layer);
  // TrkrDefs::cluskey keyhi = TrkrDefs::getClusKeyHi(trackerid, layer);
  // retpair.first = m_clusmap.lower_bound(keylo);
  // retpair.second = m_clusmap.upper_bound(keyhi);

  ConstRange retpair;
  if(phi_segment<max_phisegment&&z_segment<max_zsegment){
    retpair.first  = m_clusmap[layer][phi_segment][z_segment].begin();
    retpair.second = m_clusmap[layer][phi_segment][z_segment].end();
  }else{
    retpair.first  = m_clusmap[layer][phi_segment][z_segment].begin();
    retpair.second = m_clusmap[layer][phi_segment][z_segment].begin();
  }
  return retpair;
}

TrkrClusterContainer::Iterator
TrkrClusterContainer::findOrAddCluster(TrkrDefs::cluskey key){
  const unsigned int layer = TrkrDefs::getLayer(key);
  const unsigned int side  = TrkrDefs::getZElement(key);
  const unsigned int sector= TrkrDefs::getPhiElement(key);
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
  const unsigned int layer = TrkrDefs::getLayer(key);
  const unsigned int side  = TrkrDefs::getZElement(key);
  const unsigned int sector= TrkrDefs::getPhiElement(key);
  TrkrClusterContainer::Iterator it = m_clusmap[layer][sector][side].find(key);

  if (it != m_clusmap[layer][sector][side].end())
  {
    return it->second;
  }
  std::cout << " returning zero " << std::endl;
  return 0;
}
