#ifndef TRACKBASE_TRKRCLUSTERCROSSINGASSOCV1_H
#define TRACKBASE_TRKRCLUSTERCROSSINGASSOCV1_H
/**

 * @file trackbase/TrkrClusterCrossingAssocv1.h
 * @author Tony Frawley
 * @date March 2022
 * @brief Version 1 of class for associating clusters to the bunch crossing that created them
 */

#include "TrkrDefs.h"
#include "TrkrClusterCrossingAssoc.h"

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for pair

/**
 * @brief Class for associating clusters to the bunch crossing that created them
 *
 * Store the associations between clusters and beam crossings.
 */
class TrkrClusterCrossingAssocv1 : public TrkrClusterCrossingAssoc
{
  public:

  TrkrClusterCrossingAssocv1() = default;

  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  void addAssoc(TrkrDefs::cluskey, short int) override;

  ConstRange getCrossings(TrkrDefs::cluskey) const override;

  ConstRange getAll() const override;

  unsigned int size(void) const override;

private:

  std::multimap<TrkrDefs::cluskey, short int> m_map;

  ClassDefOverride(TrkrClusterCrossingAssocv1, 1);
};

#endif // TRACKBASE_TRKRCLUSTERCROSSINGASSOCV1_H
