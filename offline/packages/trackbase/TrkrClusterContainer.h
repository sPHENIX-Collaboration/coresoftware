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

  //! constructor
  TrkrClusterContainer() = default;

  //! reset method
  virtual void Reset() {}

  //! identify object
  virtual void identify(std::ostream &os = std::cout) const {}

  //! add a cluster
  virtual ConstIterator addCluster(TrkrCluster*);

  //! add a cluster with specific key
  virtual ConstIterator addClusterSpecifyKey(const TrkrDefs::cluskey, TrkrCluster* );

  //! remove cluster
  virtual void removeCluster(TrkrDefs::cluskey) {}

  //! remove cluster
  virtual void removeCluster(TrkrCluster* ) {}

  //! find cluster matching key if any, add a new one otherwise and return cluster
  virtual Iterator findOrAddCluster(TrkrDefs::cluskey);
  
  //! return all clusters
  virtual ConstRange getClusters() const;

  //! get all clusters matching hitset
  virtual ConstRange getClusters(TrkrDefs::hitsetkey) const;

  //! get pointer to map containing clusters mathching hitset
  virtual Map* getClusterMap(TrkrDefs::hitsetkey) { return nullptr; }

  //! find cluster matching given key
  virtual TrkrCluster* findCluster(TrkrDefs::cluskey) const { return nullptr; }

  //! total number of clusters
  virtual unsigned int size() const { return 0; }

  private:

  ClassDef(TrkrClusterContainer, 1)

};

#endif //TRACKBASE_TRKRCLUSTERCONTAINER_H
