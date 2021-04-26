/**
 * @file trackbase/TrkrClusterHitAssocv1.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Class for associating clusters to the hits that went into them
 */
#ifndef TRACKBASE_TRKRCLUSTERHITASSOCV1_H
#define TRACKBASE_TRKRCLUSTERHITASSOCV1_H

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
class TrkrClusterHitAssocv1 : public TrkrClusterHitAssoc
{
public:

  //! ctor
  TrkrClusterHitAssocv1();

  //! dtor
  virtual ~TrkrClusterHitAssocv1();

  virtual void Reset();

  virtual void identify(std::ostream &os = std::cout) const;

  /**
   * @brief Add association between cluster and hit
   * @param[in] ckey Cluster key
   * @param[in] hidx Index of the hit in TrkrHitSet
   */
  virtual void addAssoc(TrkrDefs::cluskey ckey, unsigned int hidx);

  /**
   * @brief Get all the hits associated with a cluster by key
   * @param[in] ckey Cluster key
   * @param[out] Range over hits associated with @c ckey
   */
  virtual std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator>  
    getHits(TrkrDefs::cluskey ckey);

private:

  std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey> m_map;

  ClassDef(TrkrClusterHitAssocv1, 1);
};

#endif // TRACKBASE_TRKRCLUSTERHITASSOCV1_H
