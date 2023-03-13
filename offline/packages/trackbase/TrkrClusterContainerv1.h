/**
 * @file trackbase/TrkrClusterContainerv1.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Cluster container object
 */
#ifndef TRACKBASE_TRKRCLUSTERCONTAINERV1_H
#define TRACKBASE_TRKRCLUSTERCONTAINERV1_H

#include "TrkrClusterContainer.h"
#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for pair

class TrkrCluster;

/**
 * @brief Cluster container object
 *
 * Container for TrkrCluster objects
 */
class TrkrClusterContainerv1 : public TrkrClusterContainer
{
 public:
  typedef std::map<TrkrDefs::cluskey, TrkrCluster *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  TrkrClusterContainerv1() = default;

  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  void addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster *newClus) override;

  void removeCluster(TrkrDefs::cluskey) override;

  ConstRange getClusters() const override;

  TrkrCluster *findCluster(TrkrDefs::cluskey key) const override;

  unsigned int size() const override;

 private:
  Map m_clusmap;
  ClassDefOverride(TrkrClusterContainerv1, 1)
};

#endif  // TRACKBASE_TRKRCLUSTERCONTAINER_H
