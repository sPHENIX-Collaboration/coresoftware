#ifndef TRACKBASEHISTORIC_TRACKSEEDCONTAINERV1_H
#define TRACKBASEHISTORIC_TRACKSEEDCONTAINERV1_H

#include "TrackSeed.h"
#include "TrackSeedContainer.h"

#include <iostream>
#include <vector>

class TrackSeedContainer_v1 : public TrackSeedContainer
{

 public: 
  TrackSeedContainer_v1();
  TrackSeedContainer_v1(const TrackSeedContainer_v1& trackmap);
  TrackSeedContainer_v1& operator=(const TrackSeedContainer_v1& seedContainer);
  ~TrackSeedContainer_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new TrackSeedContainer_v1(*this); }

  bool empty() const override { return m_seeds.empty(); }
  std::size_t size() const override { return m_seeds.size(); }
  void clear() override { Reset(); }

  const TrackSeed* get(const std::size_t key) const override;
  TrackSeed* get(const std::size_t key) override;
  TrackSeed* insert(const TrackSeed* seed) override;
  void erase(const std::size_t key) override
    {
      delete m_seeds.at(key);
      m_seeds.at(key) = nullptr;
    }

  std::size_t index(ConstIter it) const override { return it - m_seeds.begin(); }
  std::size_t index(Iter it) const override { return it - m_seeds.begin(); }

  ConstIter begin() const override { return m_seeds.begin(); }
  ConstIter find(const std::size_t key) const override { return m_seeds.begin() + key; }
  std::size_t find(const TrackSeed*) const override;
  ConstIter end() const override { return m_seeds.end(); }
  
  Iter begin() override { return m_seeds.begin(); }
  Iter find(const std::size_t key) override { return m_seeds.begin() + key; }
  std::size_t find(const TrackSeed*) override;
  Iter end() override { return m_seeds.end(); }

 private:
  Container m_seeds;

  ClassDefOverride(TrackSeedContainer_v1, 1);

};

#endif
