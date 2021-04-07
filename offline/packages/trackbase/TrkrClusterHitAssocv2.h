/**
 * @file trackbase/TrkrClusterHitAssocv2.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Version 1 of class for associating clusters to the hits that went into them
 */
#ifndef TRACKBASE_TRKRCLUSTERHITASSOCV2_H
#define TRACKBASE_TRKRCLUSTERHITASSOCV2_H

#include "TrkrDefs.h"
#include "TrkrClusterHitAssoc.h"

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for pair

/**
 * @brief Class for associating clusters to the hits that went into them
 *
 * Store the associations between clusters and the hits that went into them.
 */
class TrkrClusterHitAssocv2 : public TrkrClusterHitAssoc
{
public:

  //! ctor
  TrkrClusterHitAssocv2();

  //! dtor
  virtual ~TrkrClusterHitAssocv2();

  virtual void Reset();

  virtual void identify(std::ostream &os = std::cout) const;

  /**
   * @brief Add association between cluster and hit
   * @param[in] ckey Cluster key
   * @param[in] hidx Index of the hit in TrkrHitSet
   */
  virtual void addAssoc(TrkrDefs::cluskey ckey, unsigned int hidx);

  virtual std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey> *getClusterSet(unsigned int layer, unsigned int phi_segment, unsigned int z_segment);

  /**
   * @brief Get all the hits associated with a cluster by key
   * @param[in] ckey Cluster key
   * @param[out] Range over hits associated with @c ckey
   */
  virtual std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator> 
    getHits(TrkrDefs::cluskey ckey);

  virtual unsigned int size(void) const;

private:
  unsigned int max_layer = 57;
  unsigned int max_phisegment = 12;
  unsigned int max_zsegment = 8;
  std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey> m_map[57][12][8];

  ClassDef(TrkrClusterHitAssocv2, 1);
};

#endif // TRACKBASE_TRKRCLUSTERHITASSOC_H
