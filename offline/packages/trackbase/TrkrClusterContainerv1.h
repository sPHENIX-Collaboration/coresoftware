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

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

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
  
  virtual void Reset() override;

  virtual void identify(std::ostream &os = std::cout) const override;

  virtual ConstIterator addCluster(TrkrCluster *newClus) override;

  virtual ConstIterator addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster *newClus) override;

  virtual void removeCluster(TrkrDefs::cluskey) override;

  virtual void removeCluster(TrkrCluster*) override;

  virtual Iterator findOrAddCluster(TrkrDefs::cluskey key) override;
  
  virtual ConstRange getClusters() const override;

  virtual TrkrCluster *findCluster(TrkrDefs::cluskey key) const override;

  virtual unsigned int size() const override;

  private:
  Map m_clusmap;
  ClassDefOverride(TrkrClusterContainerv1, 1)
};

#endif //TRACKBASE_TRKRCLUSTERCONTAINER_H
