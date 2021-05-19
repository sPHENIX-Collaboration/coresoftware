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
  
  virtual void Reset();

  virtual void identify(std::ostream &os = std::cout) const;

  virtual ConstIterator addCluster(TrkrCluster *newClus);

  virtual ConstIterator addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster *newClus);

  virtual void removeCluster(TrkrDefs::cluskey);

  virtual void removeCluster(TrkrCluster*);

  virtual Iterator findOrAddCluster(TrkrDefs::cluskey key);
  
  virtual ConstRange getClusters() const;

  virtual TrkrCluster *findCluster(TrkrDefs::cluskey key) const;

  virtual unsigned int size() const;

  private:
  Map m_clusmap;
  ClassDef(TrkrClusterContainerv1, 1)
};

#endif //TRACKBASE_TRKRCLUSTERCONTAINER_H
