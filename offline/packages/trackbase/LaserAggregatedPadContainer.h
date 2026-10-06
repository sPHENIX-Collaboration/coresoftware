#ifndef TRACKBASE_LASERAGGREGATEDPADCONTAINER_H
#define TRACKBASE_LASERAGGREGATEDPADCONTAINER_H

/**
 * @file trackbase/LaserAggregatedPadContainer.h
 * @author Ben Kimelman
 * @date February 2024
 * @brief Laser cluster container base class
 */

#include "TrkrDefs.h"
#include "padkey.h"

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class LaserAggregatedPad;

/**
 * @brief Cluster container object
 */

class LaserAggregatedPadContainer : public PHObject
{
 public:

  //!@name convenient shortuts
  //@{
  using Map = std::map<padkey::key, LaserAggregatedPad *>;
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
  virtual void addPad(padkey::key /*pad*/, LaserAggregatedPad* /*newPad*/) = 0;

  virtual void merge(const LaserAggregatedPadContainer* /*other*/) {}

  //! remove cluster
  virtual void removePad(padkey::key /*pad*/) {}
  
  //! return all clusters
  virtual ConstRange getPads() const = 0;
  
  virtual ConstRange getPadsInBlock(uint8_t /*blk*/) const = 0;
  
  //! find cluster matching given key
  virtual LaserAggregatedPad* findPad(padkey::key /*pad*/) const { return nullptr; }

  //! total number of clusters
  virtual unsigned int size() const { return 0; }

  virtual double getNLaserEvents() const { return 0; }
  virtual void setNLaserEvents(double /*nLaser*/) {}

  protected:
  //! constructor
  LaserAggregatedPadContainer() = default;

  private:

  ClassDefOverride(LaserAggregatedPadContainer, 1)

};

#endif //TRACKBASE_LASERAGGREGATEDPADCONTAINER_H
