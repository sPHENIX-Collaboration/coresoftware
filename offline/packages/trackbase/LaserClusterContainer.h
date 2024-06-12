#ifndef TRACKBASE_LASERCLUSTERCONTAINER_H
#define TRACKBASE_LASERCLUSTERCONTAINER_H

/**
 * @file trackbase/LaserClusterContainer.h
 * @author Ben Kimelman
 * @date February 2024
 * @brief Laser cluster container base class
 */

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class LaserCluster;

/**
 * @brief Cluster container object
 */
class LaserClusterContainer : public PHObject
{
 public:

  //!@name convenient shortuts
  //@{
  using Map = std::map<TrkrDefs::cluskey, LaserCluster *>;
  using Iterator = Map::iterator;
  using ConstIterator = Map::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;
  //@}

  //! reset method
  void Reset() override {}

  //! identify object
  void identify(std::ostream &/*os*/ = std::cout) const override {}
  
  //! add a cluster with specific key
  virtual void addClusterSpecifyKey(const TrkrDefs::cluskey, LaserCluster*) = 0;

  //! remove cluster
  virtual void removeCluster(TrkrDefs::cluskey) {}
  
  //! return all clusters
  virtual ConstRange getClusters() const = 0;

  //! find cluster matching given key
  virtual LaserCluster* findCluster(TrkrDefs::cluskey) const { return nullptr; }

  //! total number of clusters
  virtual unsigned int size() const { return 0; }

  protected:
  //! constructor
  LaserClusterContainer() = default;

  private:

  ClassDefOverride(LaserClusterContainer, 1)

};

#endif //TRACKBASE_LASERCLUSTERCONTAINER_H
