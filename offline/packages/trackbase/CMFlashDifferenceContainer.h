#ifndef TRACKBASE_CMFLASHDIFFERENCECONTAINER_H
#define TRACKBASE_CMFLASHDIFFERENCECONTAINER_H

/**
 * @file trackbase/CMFlashDifferenceContainer.h
 * @author Tony Frawley
 * @date January 2022
 * @brief Central membrane flash difference container base class
 */

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class CMFlashDifference;

/**
 * @brief CM cluster difference due to distortions container object
 */
class CMFlashDifferenceContainer : public PHObject
{
 public:

  //!@name convenient shortuts
  //@{
  using Map = std::map<unsigned int, CMFlashDifference *>;
  using Iterator = Map::iterator;
  using ConstIterator = Map::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;
  //@}

  //! reset method
  void Reset() override {}

  //! identify object
  void identify(std::ostream &/*os*/ = std::cout) const override {}

  //! add a differences with specific key
  virtual void addDifferenceSpecifyKey(const unsigned int, CMFlashDifference* ) = 0;

  //! remove differences
  virtual void removeDifference(unsigned int) {}
  
  //! return all differences
  virtual ConstRange getDifferences() const = 0;

  //! find differences matching given key
  virtual CMFlashDifference* findDifference(unsigned int) const { return nullptr; }

  //! total number of differences
  virtual unsigned int size() const { return 0; }

  protected:
  //! constructor
  CMFlashDifferenceContainer() = default;

  private:

  ClassDefOverride(CMFlashDifferenceContainer, 1)

};

#endif //TRACKBASE_CMFLASHDIFFERENCECONTAINER_H
