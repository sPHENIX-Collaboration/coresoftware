#ifndef TRACKBASE_TRKRCLUSTERHITASSOCV3_H
#define TRACKBASE_TRKRCLUSTERHITASSOCV3_H
/**

 * @file trackbase/TrkrClusterHitAssocv3.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Version 3 of class for associating clusters to the hits that went into them
 */

#include "TrkrClusterHitAssoc.h"
#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for pair

/**
 * @brief Class for associating clusters to the hits that went into them
 *
 * Store the associations between clusters and the hits that went into them.
 */
class TrkrClusterHitAssocv3 : public TrkrClusterHitAssoc
{
 public:
  TrkrClusterHitAssocv3() = default;

  void Reset() override;

  void identify(std::ostream& os = std::cout) const override;

  void addAssoc(TrkrDefs::cluskey, unsigned int) override;

  Map* getClusterMap(TrkrDefs::hitsetkey) override;

  ConstRange getHits(TrkrDefs::cluskey) override;

  unsigned int size(void) const override;

 private:
  std::map<TrkrDefs::hitsetkey, Map> m_map;

  ClassDefOverride(TrkrClusterHitAssocv3, 1);
};

#endif  // TRACKBASE_TRKRCLUSTERHITASSOC_H
