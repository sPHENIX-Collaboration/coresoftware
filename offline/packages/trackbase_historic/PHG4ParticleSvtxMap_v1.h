#ifndef TRACKBASEHISTORIC_PHG4PARTICLESVTXMAP_V1_H
#define TRACKBASEHISTORIC_PHG4PARTICLESVTXMAP_V1_H

#include "PHG4ParticleSvtxMap.h"

#include <iostream>

class PHG4ParticleSvtxMap_v1 : public PHG4ParticleSvtxMap
{
 public:
  PHG4ParticleSvtxMap_v1();
  PHG4ParticleSvtxMap_v1(const PHG4ParticleSvtxMap_v1& map);
  PHG4ParticleSvtxMap_v1& operator=(const PHG4ParticleSvtxMap_v1& map);
  ~PHG4ParticleSvtxMap_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new PHG4ParticleSvtxMap_v1(*this); }
  void Reset() override { clear(); m_processed = false; }
  bool empty() const override { return m_map.empty(); }
  std::size_t size() const override { return m_map.size(); }
  std::size_t count(const int key) const override { return m_map.count(key); }
  void clear() override { m_map.clear(); }

  bool processed() const override { return m_processed; }
  void setProcessed(const bool process) override { m_processed = process; }

  const WeightedRecoTrackMap & get(const int key) const override;
  WeightedRecoTrackMap & get(const int key) override
  {
    return m_map[key];
  }
  WeightedRecoTrackMap insert(const int key, const WeightedRecoTrackMap map) override;
  std::size_t erase(const int key) override
  {
    return m_map.erase(key);
  }

  ConstIter begin() const override { return m_map.begin(); }
  ConstIter find(const int key) const override
  {
    return m_map.find(key);
  }
  ConstIter end() const override { return m_map.end(); }

  Iter begin() override { return m_map.begin(); }
  Iter find(const int key) override { return m_map.find(key); }
  Iter end() override { return m_map.end(); }

 private:
  PHG4ParticleSvtxMap::Map m_map;
  bool m_processed = false;

  ClassDefOverride(PHG4ParticleSvtxMap_v1, 1);
};

#endif
