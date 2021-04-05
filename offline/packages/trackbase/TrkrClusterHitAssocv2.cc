/**                                                                                                                                                                                                       
 * @file trackbase/TrkrClusterHitAssocv2.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief TrkrClusterHitAssoc implementation
 */

#include "TrkrClusterHitAssocv2.h"
#include "TrkrDefs.h"

#include <ostream>  // for operator<<, endl, basic_ostream, ostream, basic_o...

TrkrClusterHitAssocv2::TrkrClusterHitAssocv2() 
  : m_map()
{
}


TrkrClusterHitAssocv2::~TrkrClusterHitAssocv2()
{
}

void TrkrClusterHitAssocv2::Reset()
{
  for(unsigned int layer = 0;layer < max_layer; layer++){
    for(unsigned int phi_segment = 0;phi_segment < max_phisegment;phi_segment++){
      for(unsigned int z_segment = 0; z_segment < max_zsegment; z_segment++){
	m_map[layer][phi_segment][z_segment].clear();
	  /*
	while (m_map[layer][phi_segment][z_segment].begin() != m_map[layer][phi_segment][z_segment].end())
	  {
	    delete m_map[layer][phi_segment][z_segment].begin()->second;
	    m_map[layer][phi_segment][z_segment].erase(m_clusmap[layer][phi_segment][z_segment].begin());
	  }
	  */
      }
    }
  }
  //  m_map.clear();
  return;
}

void TrkrClusterHitAssocv2::identify(std::ostream &os) const
{
  std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator iter;
  os << "-----TrkrClusterHitAssocv2-----" << std::endl;
  os << "Number of associations: " << size() << std::endl;
  for(unsigned int layer = 0;layer < max_layer; layer++){
    for(unsigned int phi_segment = 0;phi_segment < max_phisegment;phi_segment++){
      for(unsigned int z_segment = 0; z_segment < max_zsegment; z_segment++){
	for (iter = m_map[layer][phi_segment][z_segment].begin(); 
	     iter != m_map[layer][phi_segment][z_segment].end();
	     ++iter)
	  {
	    int layer = TrkrDefs::getLayer(iter->first);
	    os << "clus key " << iter->first << std::dec
	       << " layer " << layer 
	       << " hit key: " << iter->second << std::endl;
	  }
      }
    }
  }
  /*  for ( auto& entry : m_map )
  {
    // os << "   cluster key: 0x" << std::hex << entry.first << std::dec
    os << "   cluster key: "  << entry.first << std::dec
       << " hit key: " << entry.second << std::endl;
  }
  */
  os << "------------------------------" << std::endl;

  return;

}

void 
TrkrClusterHitAssocv2::addAssoc(TrkrDefs::cluskey ckey, unsigned int hidx)
{
  unsigned int layer = TrkrDefs::getLayer(ckey);
  unsigned int sector = TrkrDefs::getPhiElement(ckey);
  unsigned int side = TrkrDefs::getZElement(ckey); 
  m_map[layer][sector][side].insert(std::make_pair(ckey, hidx));
}

std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>  
*TrkrClusterHitAssocv2::getClusterSet(unsigned int layer, unsigned int phi_segment, unsigned int z_segment)
{
  return &m_map[layer][phi_segment][z_segment];
}

std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator>
TrkrClusterHitAssocv2::getHits(TrkrDefs::cluskey ckey)
{
  unsigned int layer = TrkrDefs::getLayer(ckey);
  unsigned int sector = TrkrDefs::getPhiElement(ckey);
  unsigned int side = TrkrDefs::getZElement(ckey); 

  std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator> retpair;
  retpair.first = m_map[layer][sector][side].lower_bound(ckey);
  retpair.second = m_map[layer][sector][side].upper_bound(ckey);
  return retpair;
}

unsigned int 
TrkrClusterHitAssocv2::size(void) const
{
  unsigned int size = 0;
  for(unsigned layer = 0;layer < max_layer; layer++){
    for(unsigned phi_segment = 0;phi_segment < max_phisegment;phi_segment++){
      for(unsigned z_segment = 0; z_segment < max_zsegment; z_segment++){
	size += m_map[layer][phi_segment][z_segment].size(); 
      }
    }
  }
  return size;
}
