#ifndef TRACKBASE_CMFLASHCLUSTERCONTAINER_H
#define TRACKBASE_CMFLASHCLUSTERCONTAINER_H

/**
 * @file trackbase/CMFlashClusterContainer.h
 * @author Tony Frawley
 * @date January 2022
 * @brief Central membrane flash cluster container base class
 */

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class CMFlashCluster;

/**
 * @brief Cluster container object
 */
class CMFlashClusterContainer : public PHObject
{
 public:

  //!@name convenient shortuts
  //@{
  using Map = std::map<unsigned int, CMFlashCluster *>;
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
  virtual void addClusterSpecifyKey(const unsigned int, CMFlashCluster* ) = 0;

  //! remove cluster
  virtual void removeCluster(unsigned int) {}
  
  //! return all clusters
  virtual ConstRange getClusters() const = 0;

  //! find cluster matching given key
  virtual CMFlashCluster* findCluster(unsigned int) const { return nullptr; }

  //! total number of clusters
  virtual unsigned int size() const { return 0; }

  protected:
  //! constructor
  CMFlashClusterContainer() = default;

  private:

  ClassDefOverride(CMFlashClusterContainer, 1)

};

#endif //TRACKBASE_CMFLASHCLUSTERCONTAINER_H
