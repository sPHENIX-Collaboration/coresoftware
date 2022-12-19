#ifndef TRACKBASEHISTORIC_SVTXALIGNMENTSTATEMAP_V1_H
#define TRACKBASEHISTORIC_SVTXALIGNMENTSTATEMAP_V1_H

#include "SvtxAlignmentState.h"
#include "SvtxAlignmentStateMap.h"

class PHObject;

class SvtxAlignmentStateMap_v1 : public SvtxAlignmentStateMap
{
 public:
  SvtxAlignmentStateMap_v1();
  ~SvtxAlignmentStateMap_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new SvtxAlignmentStateMap_v1(*this); }

  bool empty() const override { return m_map.empty(); }
  std::size_t size() const override { return m_map.size(); }
  std::size_t count(unsigned int track) const override { return m_map.count(track); }
  void clear() override { Reset(); }

  const StateVec get(unsigned int track) const override;
  StateVec get(unsigned int track) override;
  StateVec insertWithKey(unsigned int track, StateVec states) override;
  std::size_t erase(unsigned int track) override;

  ConstIter begin() const override { return m_map.begin(); }
  ConstIter find(unsigned int track) const override { return m_map.find(track); }
  ConstIter end() const override { return m_map.end(); }

  Iter begin() override { return m_map.begin(); }
  Iter find(unsigned int track) override { return m_map.find(track); }
  Iter end() override { return m_map.end(); }

 private:
  StateMap m_map;

  ClassDefOverride(SvtxAlignmentStateMap_v1, 1);
};

#endif
