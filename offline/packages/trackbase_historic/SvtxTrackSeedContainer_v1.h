#ifndef TRACKBASEHISTORIC_SVTXTRACKSEEDCONTAINERV1_H
#define TRACKBASEHISTORIC_SVTXTRACKSEEDCONTAINERV1_H

#include "SvtxTrackSeed.h"
#include "SvtxTrackSeedContainer.h"

#include <iostream>
#include <vector>

class SvtxTrackSeedContainer_v1 : public SvtxTrackSeedContainer
{

 public: 
  SvtxTrackSeedContainer_v1();
  SvtxTrackSeedContainer_v1(const SvtxTrackSeedContainer_v1& trackmap);
  SvtxTrackSeedContainer_v1& operator=(const SvtxTrackSeedContainer_v1& seedContainer);
  ~SvtxTrackSeedContainer_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new SvtxTrackSeedContainer_v1(*this); }

  bool empty() const override { return m_seeds.empty(); }
  std::size_t size() const override { return m_seeds.size(); }
  void clear() override { Reset(); }

  const SvtxTrackSeed* get(const std::size_t key) const override;
  SvtxTrackSeed* get(const std::size_t key) override;
  SvtxTrackSeed* insert(const SvtxTrackSeed* seed) override;
  Iter erase(const std::size_t key) override
    {
      delete m_seeds.at(key);
      Iter iter = m_seeds.begin() + key;
      return m_seeds.erase(iter);
    }

  ConstIter begin() const override { return m_seeds.begin(); }
  ConstIter find(const std::size_t key) const override { return m_seeds.begin() + key; }
  ConstIter end() const override { return m_seeds.end(); }
  
  Iter begin() override { return m_seeds.begin(); }
  Iter find(const std::size_t key) override { return m_seeds.begin() + key; }
  Iter end() override { return m_seeds.end(); }

 private:
  TrackSeedContainer m_seeds;

  ClassDefOverride(SvtxTrackSeedContainer_v1, 1);

};

#endif
