#ifndef TRACKBASE_TRKRCLUSTERITERATIONV1_H
#define TRACKBASE_TRKRCLUSTERITERATIONV1_H
/**

 * @file trackbase/TrkrClusterHitAssocv3.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Version 3 of class for associating clusters to the hits that went into them
 */

#include "TrkrClusterIterationMap.h"
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
class TrkrClusterIterationMapv1 : public TrkrClusterIterationMap
{
 public:
  TrkrClusterIterationMapv1() = default;

  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  void addIteration(TrkrDefs::cluskey, short int) override;

  short int getIteration(TrkrDefs::cluskey ckey) override;

  unsigned int size(void) const override;

  ConstIter begin() const override { return m_map.begin(); }
  ConstIter end() const override { return m_map.end(); }

 private:
  Map m_map;

  ClassDefOverride(TrkrClusterIterationMapv1, 1);
};

#endif  // TRACKBASE_TRKRCLUSTERITERATIONV1_H
