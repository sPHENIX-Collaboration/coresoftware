#ifndef TRACKBASEHISTORIC_SVTXALIGNMENTSTATEMAP_V1_H
#define TRACKBASEHISTORIC_SVTXALIGNMENTSTATEMAP_V1_H

#include "SvtxAlignmentState.h"
#include "SvtxAlignmentStateMap.h"

class PHObject;

class SvtxAlignmentStateMap_v1 : public SvtxAlignmentStateMap
{
 public:
  SvtxAlignmentStateMap_v1();
  SvtxAlignmentStateMap_v1(const SvtxAlignmentStateMap_v1& map);
  SvtxAlignmentStateMap_v1& operator=(const SvtxAlignmentStateMap_v1& map);
  ~SvtxAlignmentStateMap_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new SvtxAlignmentStateMap_v1(*this); }
  
  bool empty() const override { return m_map.empty(); }
  std::size_t size() const override { return m_map.size(); }
  std::size_t count(TrkrDefs::cluskey ckey) const override { return m_map.count(ckey); }
  void clear() override { Reset(); }

  const SvtxAlignmentState* get(TrkrDefs::cluskey ckey) const override;
  SvtxAlignmentState* get(TrkrDefs::cluskey ckey) override;
  SvtxAlignmentState* insertWithKey(TrkrDefs::cluskey ckey, SvtxAlignmentState* state) override;
  virtual std::size_t erase(TrkrDefs::cluskey ckey) override
  {
    delete m_map[ckey];
    return m_map.erase(ckey);
  }

  ConstIter begin() const override { return m_map.begin(); }
  ConstIter find(TrkrDefs::cluskey ckey) const override { return m_map.find(ckey); }
  ConstIter end() const override { return m_map.end(); }

  Iter begin() override { return m_map.begin(); }
  Iter find(TrkrDefs::cluskey ckey) override { return m_map.find(ckey); }
  Iter end() override { return m_map.end(); }

 private:
  StateMap m_map;

  ClassDefOverride(SvtxAlignmentStateMap_v1, 1);
};


#endif
