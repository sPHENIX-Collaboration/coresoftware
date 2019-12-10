/**
 * @file trackbase/TrkrHitCellAssoc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Association object for PHG4Cells contributiong to TrkrHits
 */
#ifndef TRACKBASE_TRKRHITCELLASSOC_H
#define TRACKBASE_TRKRHITCELLASSOC_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <g4detectors/PHG4CellDefs.h>

#include <iostream>                    // for cout, ostream
#include <map>
#include <utility>                     // for pair

/**
 * @brief Association object for PHG4Cells contributiong to TrkrHits
 *
 * Association object holding a multimap of PHG4Cells associated with a given TrkrHit
 */
class TrkrHitCellAssoc : public PHObject
{
public:
  //! typedefs for convenience                                                                                                                                                                             
  typedef std::multimap<TrkrDefs::cluskey, PHG4CellDefs::keytype> MMap;
  typedef MMap::iterator Iterator;
  typedef MMap::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  //! ctor                                                                                                                                                                                                  
  TrkrHitCellAssoc();
  //! dtor                                                                                                                                                                                                  
  virtual ~TrkrHitCellAssoc();

  void Reset();

  void identify(std::ostream &os = std::cout) const;

  /**
   * @brief Add an association between hit and cell
   * @param[in] hset TrkrHitSet key
   * @param[in] hidx TrkrHit index in TrkrHitSet
   * @param[in] ckey Key for assocuated PHG4Cell
   */
  void addAssoc(const TrkrDefs::hitsetkey hset, const unsigned int hidx, const PHG4CellDefs::keytype ckey);

  /**
   * @brief Get cell keys associated with desired hit
   * @param[in] hset TrkrHitSet key
   * @param[in] hidx TrkrHit index in TrkrHitSet
   */
  ConstRange getCells(const TrkrDefs::hitsetkey hset, const unsigned int hidx);

private:
  MMap m_map;
  ClassDef(TrkrHitCellAssoc, 1);

};

#endif //TRACKBASE_TRKRHITCELLASSOC_H
