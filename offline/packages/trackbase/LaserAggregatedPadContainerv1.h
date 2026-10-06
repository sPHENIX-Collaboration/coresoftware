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
  using Map = std::map<padkey::key, LaserAggregatedPad *>;
  using Iterator = Map::iterator;
  using ConstIterator = Map::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  LaserAggregatedPadContainerv1() = default;
  ~LaserAggregatedPadContainerv1() override { LaserAggregatedPadContainerv1::Reset(); }
  LaserAggregatedPadContainerv1(const LaserAggregatedPadContainerv1&) = delete;
  LaserAggregatedPadContainerv1& operator=(const LaserAggregatedPadContainerv1&) = delete;
  
  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  void addPad(padkey::key pad, LaserAggregatedPad *newPad) override;

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
