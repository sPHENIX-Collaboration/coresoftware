#ifndef TRACKBASE_TRKRHITTRUTHASSOC_H
#define TRACKBASE_TRKRHITTRUTHASSOC_H

/**
 * @file trackbase/TrkrHitCellAssoc
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * @brief Association object for PHG4Hits contributiong to TrkrHits
 */

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
  using MMap = std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >; 
  using Iterator = MMap::iterator;
  using ConstIterator = MMap::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  void Reset() override
  {}

  void identify(std::ostream &/*os*/ = std::cout) const override
  {}

  /**
   * @brief Add an association between hit and cell
   * @param[in] hset TrkrHitSet key
   * @param[in] hidx TrkrHit index in TrkrHitSet
   * @param[in] ckey Key for assocuated g4hit
   */
  virtual void addAssoc(const TrkrDefs::hitsetkey /*hitsetkey*/, const TrkrDefs::hitkey /*hitkey*/, const PHG4HitDefs::keytype  /*g4hitkey*/) 
  {}

  /**
   * @brief Add an association between hit and cell if it does not already exist
   * @param[in] hset TrkrHitSet key
   * @param[in] hidx TrkrHit index in TrkrHitSet
   * @param[in] ckey Key for assocuated g4hit
   */
  virtual void findOrAddAssoc(const TrkrDefs::hitsetkey /*hitsetkey*/, const TrkrDefs::hitkey /*hitkey*/, const PHG4HitDefs::keytype  /*g4hitkey*/)
  {}

  virtual void removeAssoc(const TrkrDefs::hitsetkey /*hitsetkey*/, const TrkrDefs::hitkey /*hitkey*/)
  {}

  /**
   * @brief Get cell keys associated with desired hit
   * @param[in] hset TrkrHitSet key
   * @param[in] hidx TrkrHit index in TrkrHitSet
   */
  virtual void getG4Hits(const TrkrDefs::hitsetkey /*hitsetkey*/, const unsigned int /*hidx*/, MMap &/*temp_map*/) const
  {}

  protected:
  //! ctor
  TrkrHitTruthAssoc() = default;

  private:

  ClassDefOverride(TrkrHitTruthAssoc, 1);

};

#endif //TRACKBASE_TRKRHITTRUTHASSOC_H
