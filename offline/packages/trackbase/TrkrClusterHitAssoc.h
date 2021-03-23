/**
 * @file trackbase/TrkrClusterHitAssoc.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Class for associating clusters to the hits that went into them
 */
#ifndef TRACKBASE_TRKRCLUSTERHITASSOC_H
#define TRACKBASE_TRKRCLUSTERHITASSOC_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for pair

/**
 * @brief Class for associating clusters to the hits that went into them
 *
 * Store the associations between clusters and the hits that went into them.
 */
class TrkrClusterHitAssoc : public PHObject
{
public:
  //! typedefs for convenience
  typedef std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey> MMap;
  typedef MMap::iterator Iterator;
  typedef MMap::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  //! ctor
  TrkrClusterHitAssoc();
  //! dtor
  virtual ~TrkrClusterHitAssoc();

  void Reset();

  void identify(std::ostream &os = std::cout) const;

  /**
   * @brief Add association between cluster and hit
   * @param[in] ckey Cluster key
   * @param[in] hidx Index of the hit in TrkrHitSet
   */
  void addAssoc(TrkrDefs::cluskey ckey, unsigned int hidx);
  std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey> *getClusterSet(unsigned int layer, unsigned int phi_segment, unsigned int z_segment){
    return &m_map[layer][phi_segment][z_segment];
  }
  /**
   * @brief Get all the hits associated with a cluster by key
   * @param[in] ckey Cluster key
   * @param[out] Range over hits associated with @c ckey
   */
  ConstRange getHits(TrkrDefs::cluskey ckey);

  unsigned int size(void) const
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

private:
  unsigned int max_layer = 57;
  unsigned int max_phisegment = 12;
  unsigned int max_zsegment = 8;
  MMap m_map[57][12][8];
  ClassDef(TrkrClusterHitAssoc, 1);
};

#endif // TRACKBASE_TRKRCLUSTERHITASSOC_H
