#ifndef TRACKBASE_TRKRCLUSTERCONTAINER_H
#define TRACKBASE_TRKRCLUSTERCONTAINER_H

/**
 * @file trackbase/TrkrClusterContainer.h
 * @author D. McGlinchey, Hugo Pereira Da Costa
 * @date June 2018
 * @brief Cluster container base class
 */

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class TrkrCluster;

/**
 * @brief Cluster container object
 */
class TrkrClusterContainer : public PHObject
{
 public:

  //!@name convenient shortuts
  //@{
  using Map = std::map<TrkrDefs::cluskey, TrkrCluster *>;
  using Iterator = Map::iterator;
  using ConstIterator = Map::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;
  //@}

  //! reset method
  void Reset() override {}

  //! identify object
  void identify(std::ostream &/*os*/ = std::cout) const override {}

  //! add a cluster
  virtual void addCluster(TrkrCluster*) {}

  //! add a cluster with specific key
  virtual void addClusterSpecifyKey(const TrkrDefs::cluskey, TrkrCluster* ) {}

  //! remove cluster
  virtual void removeCluster(TrkrDefs::cluskey) {}

  //! remove cluster
  virtual void removeCluster(TrkrCluster* ) {}
  
  //! return all clusters
  virtual ConstRange getClusters() const;

  //! get all clusters matching hitset
  virtual ConstRange getClusters(TrkrDefs::hitsetkey) const;

  //! find cluster matching given key
  virtual TrkrCluster* findCluster(TrkrDefs::cluskey) const { return nullptr; }

  //! total number of clusters
  virtual unsigned int size() const { return 0; }

  protected:
  //! constructor
  TrkrClusterContainer() = default;

  private:

  ClassDefOverride(TrkrClusterContainer, 1)

};

#endif //TRACKBASE_TRKRCLUSTERCONTAINER_H
