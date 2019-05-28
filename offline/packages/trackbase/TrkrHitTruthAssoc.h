/**
 * @file trackbase/TrkrHitCellAssoc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Association object for PHG4Hits contributiong to TrkrHits
 */
#ifndef TRACKBASE_TRKRHITTRUTHASSOC_H
#define TRACKBASE_TRKRHITTRUTHASSOC_H

#include "TrkrDefs.h"

#include <g4main/PHG4HitDefs.h>

#include <phool/PHObject.h>

#include <iostream>              // for cout, ostream
#include <map>
#include <utility>               // for pair

/**
 * @brief Association object for PHG4Cells contributiong to TrkrHits
 *
 * Association object holding a multimap of PHG4Cells associated with a given TrkrHit
 */
class TrkrHitTruthAssoc : public PHObject
{
public:
  //! typedefs for convenience 
  typedef std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > MMap; 
  typedef MMap::iterator Iterator;
  typedef MMap::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  //! ctor                                                                                                                                                                                                  
  TrkrHitTruthAssoc();
  //! dtor                                                                                                                                                                                                  
  virtual ~TrkrHitTruthAssoc();

  void Reset();

  void identify(std::ostream &os = std::cout) const;

  /**
   * @brief Add an association between hit and cell
   * @param[in] hset TrkrHitSet key
   * @param[in] hidx TrkrHit index in TrkrHitSet
   * @param[in] ckey Key for assocuated g4hit
   */
  void addAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey, const PHG4HitDefs::keytype  g4hitkey);

  /**
   * @brief Add an association between hit and cell if it does not already exist
   * @param[in] hset TrkrHitSet key
   * @param[in] hidx TrkrHit index in TrkrHitSet
   * @param[in] ckey Key for assocuated g4hit
   */
  void findOrAddAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey, const PHG4HitDefs::keytype  g4hitkey);

  void removeAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey);

  /**
   * @brief Get cell keys associated with desired hit
   * @param[in] hset TrkrHitSet key
   * @param[in] hidx TrkrHit index in TrkrHitSet
   */
  void  getG4Hits(const TrkrDefs::hitsetkey hitsetkey, const unsigned int hidx, MMap &temp_map);

private:
  MMap m_map;
  ClassDef(TrkrHitTruthAssoc, 1);

};

#endif //TRACKBASE_TRKRHITTRUTHASSOC_H
