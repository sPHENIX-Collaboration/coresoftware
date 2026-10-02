/**
 * @file trackbase/LaserAggregatedPadContainerv1.h
 * @author Ben Kimelman
 * @date February 2024
 * @brief Implementation of laser cluster container object
 */
#ifndef TRACKBASE_LASERAGGREGATEDPADCONTAINERV1_H
#define TRACKBASE_LASERAGGREGATEDPADCONTAINERV1_H

#include "LaserAggregatedPadContainer.h"

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class LaserAggregatedPad;

/**
 * @brief laser cluster container object
 *
 * Container for LaserCluster objects
 */
class LaserAggregatedPadContainerv1 : public LaserAggregatedPadContainer
{
 public:
  typedef std::map<padkey::key, LaserAggregatedPad *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  LaserAggregatedPadContainerv1() = default;
  
  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  void addPad(const padkey::key pad, LaserAggregatedPad *newPad) override;

  void merge(const LaserAggregatedPadContainer *other) override;

  void removePad(padkey::key pad) override;
  
  ConstRange getPads() const override;

  ConstRange getPadsInBlock(uint8_t blk) const override
  {
    return std::make_pair(m_padMap.lower_bound(padkey::block_begin(blk)),
                        m_padMap.lower_bound(padkey::block_end(blk)));
  }

  LaserAggregatedPad *findPad(padkey::key pad) const override;

  unsigned int size() const override;

  double getNLaserEvents() const override { return nLaserEvents; }
  void setNLaserEvents(double nLaser) override { nLaserEvents = nLaser; }



  private:
  Map m_padMap;
  double nLaserEvents {0.0};
  ClassDefOverride(LaserAggregatedPadContainerv1, 1)
};

#endif //TRACKBASE_LASERAGGREGATEDPADCONTAINERV1_H
