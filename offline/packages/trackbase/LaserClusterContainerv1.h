/**
 * @file trackbase/LaserClusterContainerv1.h
 * @author Ben Kimelman
 * @date February 2024
 * @brief Implementation of laser cluster container object
 */
#ifndef TRACKBASE_LASERCLUSTERCONTAINERV1_H
#define TRACKBASE_LASERCLUSTERCONTAINERV1_H

#include "LaserClusterContainer.h"

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class LaserCluster;

/**
 * @brief laser cluster container object
 *
 * Container for LaserCluster objects
 */
class LaserClusterContainerv1 : public LaserClusterContainer
{
 public:
  typedef std::map<unsigned int, LaserCluster *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  LaserClusterContainerv1() = default;
  
  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  void addClusterSpecifyKey(const unsigned int, LaserCluster *newClus) override;

  void removeCluster(unsigned int) override;
  
  ConstRange getClusters() const override;

  LaserCluster *findCluster(unsigned int key) const override;

  unsigned int size() const override;

  private:
  Map m_clusmap;
  ClassDefOverride(LaserClusterContainerv1, 1)
};

#endif //TRACKBASE_LASERCLUSTERCONTAINER_H
