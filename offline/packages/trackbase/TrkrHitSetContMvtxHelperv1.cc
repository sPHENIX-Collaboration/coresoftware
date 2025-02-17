/**
 * @file trackbase/TrkrHitSetContMvtxHelper.cc
 * @author Yasser Corrales Morales <ycmorales@bnl.gov>
 * @date Febraury 2025
 * base class for Mvtx hitsetkey container per strobe
 */

#include "TrkrHitSetContMvtxHelperv1.h"
#include "TrkrHitSetContMvtxHelper.h"

void TrkrHitSetContMvtxHelperv1::Reset()
{
  for (auto&& [strobe, hitsetkey_set] : m_strb_hitsetkey_map)
  {
    hitsetkey_set.clear();
  }
  m_strb_hitsetkey_map.clear();
}

bool TrkrHitSetContMvtxHelperv1::addHitSetKey(const int32_t& strobe, const TrkrDefs::hitsetkey& key)
{
  const auto ret = m_strb_hitsetkey_map[strobe].insert(key);
  return ret.second;
}

const TrkrHitSetContMvtxHelper::map_tp_value&
TrkrHitSetContMvtxHelperv1::getHitSetKeys(const int32_t strobe)
{
  return m_strb_hitsetkey_map[strobe];
}
