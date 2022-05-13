#ifndef TRACKBASEHISTORIC_SVTXPHG4PARTICLEMAP_V1_H
#define TRACKBASEHISTORIC_SVTXPHG4PARTICLEMAP_V1_H

#include "SvtxPHG4ParticleMap.h"
#include <iostream>

class SvtxPHG4ParticleMap_v1 : public SvtxPHG4ParticleMap
{
 public:
  SvtxPHG4ParticleMap_v1();
  SvtxPHG4ParticleMap_v1(const SvtxPHG4ParticleMap_v1& map);
  SvtxPHG4ParticleMap_v1& operator=(const SvtxPHG4ParticleMap_v1& map);
  ~SvtxPHG4ParticleMap_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new SvtxPHG4ParticleMap_v1(*this); }
  void Reset() override { clear(); m_processed = false; }

  bool empty() const override { return m_map.empty(); }
  std::size_t size() const override { return m_map.size(); }
  std::size_t count(const unsigned int key) const override { return m_map.count(key); }
  void clear() override { m_map.clear(); }

  bool processed() const override { return m_processed; }
  void setProcessed(const bool process) override { m_processed = process; }

  const WeightedTruthTrackMap & get(const unsigned int key) const override;
  WeightedTruthTrackMap & get(const unsigned int key) override
  {
    return m_map[key];
  }
  WeightedTruthTrackMap insert(const unsigned int key, const WeightedTruthTrackMap map) override;
  std::size_t erase(const unsigned int key) override
  {
    return m_map.erase(key);
  }

  ConstIter begin() const override { return m_map.begin(); }
  ConstIter find(const unsigned int key) const override
  {
    return m_map.find(key);
  }
  ConstIter end() const override { return m_map.end(); }

  Iter begin() override { return m_map.begin(); }
  Iter find(const unsigned int key) override { return m_map.find(key); }
  Iter end() override { return m_map.end(); }

 private:
  SvtxPHG4ParticleMap::Map m_map;
  bool m_processed = false;

  ClassDefOverride(SvtxPHG4ParticleMap_v1, 1);
};

#endif
