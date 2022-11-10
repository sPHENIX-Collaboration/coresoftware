/**
 * @file trackbase/TrkrClusterHitAssocv2.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Version 1 of class for associating clusters to the hits that went into them
 */
#ifndef TRACKBASE_TRKRCLUSTERHITASSOCV2_H
#define TRACKBASE_TRKRCLUSTERHITASSOCV2_H

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
class TrkrClusterHitAssocv2 : public TrkrClusterHitAssoc
{
 public:
  TrkrClusterHitAssocv2() = default;

  void Reset() override;

  void identify(std::ostream& os = std::cout) const override;

  void addAssoc(TrkrDefs::cluskey, unsigned int) override;

  Map* getClusterMap(TrkrDefs::hitsetkey) override;

  ConstRange getHits(TrkrDefs::cluskey) override;

  unsigned int size() const override;

 private:
  unsigned int max_layer = 57;
  unsigned int max_phisegment = 12;
  unsigned int max_zsegment = 8;
  Map m_map[57][12][8];

  ClassDefOverride(TrkrClusterHitAssocv2, 1);
};

#endif  // TRACKBASE_TRKRCLUSTERHITASSOC_H
