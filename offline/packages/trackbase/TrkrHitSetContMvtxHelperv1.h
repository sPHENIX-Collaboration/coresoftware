/**
 * @file trackbase/TrkrHitSetContMvtxHelperv1.h
 * @author Yasser Corrales Morales <ycmorales@bnl.gov>
 * @date Febraury 2025
 * base class for Mvtx hitsetkey container per strobe
 */

#ifndef TRACKBASE_TRKRHITSETCONTMVTXHELPERV1_H
#define TRACKBASE_TRKRHITSETCONTMVTXHELPERV1_H

#include "TrkrHitSetContMvtxHelper.h"

/**
 * Container for Mvtx hitsekey per strobe
 */
class TrkrHitSetContMvtxHelperv1 : public TrkrHitSetContMvtxHelper
{
 public:
  //! ctor
  TrkrHitSetContMvtxHelperv1() = default;
  //! cp/mv ctor
  TrkrHitSetContMvtxHelperv1(const TrkrHitSetContMvtxHelperv1 &) = default;
  TrkrHitSetContMvtxHelperv1(TrkrHitSetContMvtxHelperv1 &&) = default;
  //! cp/mv assignment
  TrkrHitSetContMvtxHelperv1 &operator=(const TrkrHitSetContMvtxHelperv1 &) = default;
  TrkrHitSetContMvtxHelperv1 &operator=(TrkrHitSetContMvtxHelperv1 &&) = default;

  //! dtor
  ~TrkrHitSetContMvtxHelperv1() override = default;

  //! PHObject functions
  void Reset() override;

  bool addHitSetKey(const int32_t &strobe, const TrkrDefs::hitsetkey &hitsetkey) override;

  //! preferred removal method, key is currently the hit id
  bool removeHitSetKey(const int32_t strobe, const TrkrDefs::hitsetkey hitsetkey) override
  {
    size_t ret = m_strb_hitsetkey_map[strobe].erase(hitsetkey);
    if (!ret)
    {
      std::cout << "hitsetkey " << hitsetkey << " was not found in strobe " << strobe << std::endl;
      return false;
    }
    return true;
  }

  //! return all HitSetKeys
  const TrkrHitSetContMvtxHelper::map_tp_value &getHitSetKeys(int32_t strobe) override;

  unsigned int size() const override
  {
    return 0;
  }

 protected:
 private:
  Map m_strb_hitsetkey_map;

  ClassDefOverride(TrkrHitSetContMvtxHelperv1, 1)
};

#endif  // TRACKBASE_TRKRHITSETCONTAINER_H
