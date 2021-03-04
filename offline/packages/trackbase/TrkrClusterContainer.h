/**
 * @file trackbase/TrkrClusterContainer.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Cluster container object
 */
#ifndef TRACKBASE_TRKRCLUSTERCONTAINER_H
#define TRACKBASE_TRKRCLUSTERCONTAINER_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class TrkrCluster;

/**
 * @brief Cluster container object
 *
 * Container for TrkrCluster objects
 */
class TrkrClusterContainer : public PHObject
{
 public:
  typedef std::map<TrkrDefs::cluskey, TrkrCluster *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  TrkrClusterContainer(){}

  virtual ~TrkrClusterContainer() {}
  void Reset();

  void identify(std::ostream &os = std::cout) const;

  ConstIterator addCluster(TrkrCluster *newClus);
  ConstIterator addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster *newClus);

  //! preferred removal method, key is currently the clus id
  void removeCluster(TrkrDefs::cluskey key);

  //! inefficent, use key where possible instead
  void removeCluster(TrkrCluster *clus);

  Iterator findOrAddCluster(TrkrDefs::cluskey key);

  //! return all Clusters matching a given detid
  ConstRange getClusters(const TrkrDefs::TrkrId trackerid) const;

  //! return all Clusters matching a given detid and layer
  ConstRange getClusters(const TrkrDefs::TrkrId trackerid, const unsigned int layer) const;

  //! return all Clusters in a given layer
  ConstRange getClusters(const unsigned int layer) const;

  //! return all clusters
  ConstRange getClusters(void) const;

  //! return all clusters
  std::map<TrkrDefs::cluskey, TrkrCluster *> *getClusterSet(unsigned int layer, unsigned int phi_segment, unsigned int z_segment){
    return &m_clusmap[layer][phi_segment][z_segment];
  }

  TrkrCluster *findCluster(TrkrDefs::cluskey key);

  unsigned int size(void) const
  {
    unsigned int size = 0;
    for(unsigned layer = 0;layer < max_layer; layer++){
      for(unsigned phi_segment = 0;phi_segment < max_phisegment;phi_segment++){
	for(unsigned z_segment = 0; z_segment < max_zsegment; z_segment++){
	  size += m_clusmap[layer][phi_segment][z_segment].size(); 
	}
      }
    }
    return size;
  }

 protected:
  unsigned int max_layer = 57;
  unsigned int max_phisegment = 12;
  unsigned int max_zsegment = 8;
  Map m_clusmap[57][12][8];
  ClassDef(TrkrClusterContainer, 1)
};

#endif //TRACKBASE_TRKRCLUSTERCONTAINER_H
