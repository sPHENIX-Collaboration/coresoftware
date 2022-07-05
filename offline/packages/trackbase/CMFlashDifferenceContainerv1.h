/**
 * @file trackbase/CMFlashDifferenceContainerv1.h
 * @author Tony Frawley
 * @date January 2022
 * @brief Implementation of central membrane flash difference container object
 */
#ifndef TRACKBASE_CMFLASHDIFFERENCECONTAINERV1_H
#define TRACKBASE_CMFLASHDIFFERENCECONTAINERV1_H

#include "CMFlashDifferenceContainer.h"

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class CMFlashDifference;

/**
 * @brief CM flash difference container object
 *
 * Container for CMFlashDifference objects
 */
class CMFlashDifferenceContainerv1 : public CMFlashDifferenceContainer
{
 public:
  typedef std::map<unsigned int, CMFlashDifference *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  CMFlashDifferenceContainerv1() = default;
  
  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  void addDifferenceSpecifyKey(const unsigned int, CMFlashDifference *newClus) override;

  void removeDifference(unsigned int) override;
  
  ConstRange getDifferences() const override;

  CMFlashDifference *findDifference(unsigned int key) const override;

  unsigned int size() const override;

  private:
  Map m_map;
  ClassDefOverride(CMFlashDifferenceContainerv1, 1)
};

#endif //TRACKBASE_CMFLASHDIFFERENCECONTAINER_H
