/**
 * @file trackbase/CMFlashClusterContainerv1.h
 * @author Tony Frawley
 * @date January 2022
 * @brief Implementation of central membrane flash cluster container object
 */
#ifndef TRACKBASE_CMFLASHCLUSTERCONTAINERV1_H
#define TRACKBASE_CMFLASHCLUSTERCONTAINERV1_H

#include "CMFlashClusterContainer.h"

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class CMFlashCluster;

/**
 * @brief CM flash cluster container object
 *
 * Container for CMFlashCluster objects
 */
class CMFlashClusterContainerv1 : public CMFlashClusterContainer
{
 public:
  typedef std::map<unsigned int, CMFlashCluster *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  CMFlashClusterContainerv1() = default;
  
  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  void addClusterSpecifyKey(const unsigned int, CMFlashCluster *newClus) override;

  void removeCluster(unsigned int) override;
  
  ConstRange getClusters() const override;

  CMFlashCluster *findCluster(unsigned int key) const override;

  unsigned int size() const override;

  private:
  Map m_clusmap;
  ClassDefOverride(CMFlashClusterContainerv1, 1)
};

#endif //TRACKBASE_CMFLASHCLUSTERCONTAINER_H
