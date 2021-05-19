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

  TrkrClusterHitAssocv2() = default;

  virtual void Reset();

  virtual void identify(std::ostream &os = std::cout) const;

  virtual void addAssoc(TrkrDefs::cluskey, unsigned int);

  virtual Map* getClusterMap(TrkrDefs::hitsetkey);

  virtual ConstRange getHits(TrkrDefs::cluskey);

  virtual unsigned int size() const;

private:
  unsigned int max_layer = 57;
  unsigned int max_phisegment = 12;
  unsigned int max_zsegment = 8;
  Map m_map[57][12][8];

  ClassDef(TrkrClusterHitAssocv2, 1);
};

#endif // TRACKBASE_TRKRCLUSTERHITASSOC_H
