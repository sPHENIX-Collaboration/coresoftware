#ifndef TRACKBASE_TRKRHITTRUTHASSOCV1_H
#define TRACKBASE_TRKRHITTRUTHASSOCV1_H

/**
 * @file trackbase/TrkrHitCellAssoc
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * @brief Association object for PHG4Hits contributiong to TrkrHits
 */

#include "TrkrDefs.h"
#include "TrkrHitTruthAssoc.h"

#include <iostream>              // for cout, ostream
#include <map>
#include <utility>               // for pair

/**
 * @brief Association object for PHG4Cells contributiong to TrkrHits
 *
 * Association object holding a multimap of PHG4Cells associated with a given TrkrHit
 */
class TrkrHitTruthAssocv1 : public TrkrHitTruthAssoc
{
  
  public:

  TrkrHitTruthAssocv1() = default;
 
  virtual void Reset();

  virtual void identify(std::ostream &os = std::cout) const;

  virtual void addAssoc(const TrkrDefs::hitsetkey, const TrkrDefs::hitkey, const PHG4HitDefs::keytype);

  virtual void findOrAddAssoc(const TrkrDefs::hitsetkey, const TrkrDefs::hitkey, const PHG4HitDefs::keytype);

  virtual void removeAssoc(const TrkrDefs::hitsetkey, const TrkrDefs::hitkey);
  
  virtual void getG4Hits(const TrkrDefs::hitsetkey hitsetkey, const unsigned int hidx, MMap &temp_map) const;

  private:
  
  MMap m_map;
  
  ClassDef(TrkrHitTruthAssocv1, 1);

};

#endif //TRACKBASE_TRKRHITTRUTHASSOC_H
