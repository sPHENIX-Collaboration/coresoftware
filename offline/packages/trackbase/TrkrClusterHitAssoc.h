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

#include <map>

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

  /**
   * @brief Get all the hits associated with a cluster by key
   * @param[in] ckey Cluster key
   * @param[out] Range over hits associated with @c ckey
   */
  ConstRange getHits(TrkrDefs::cluskey ckey);

private:
  MMap m_map;
  ClassDef(TrkrClusterHitAssoc, 1);
};

#endif // TRACKBASE_TRKRCLUSTERHITASSOC_H
