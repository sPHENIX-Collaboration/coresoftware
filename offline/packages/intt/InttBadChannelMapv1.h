#ifndef INTT_BAD_CHANNEL_MAPv1_H
#define INTT_BAD_CHANNEL_MAPv1_H

#include "InttMap.h"
#include "InttBadChannelMap.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <set>

class CDBTTree;

class InttBadChannelMapv1 : public InttBadChannelMap
{
 public:
  InttBadChannelMapv1() = default;
  ~InttBadChannelMapv1() override;

  void identify(std::ostream& = std::cout) const override;
  std::size_t size() const override;

  bool IsBad(InttMap::Offline_s const&) const override;

 protected:
  int v_LoadFromCDBTTree(CDBTTree&) override;

 private:
  typedef std::set<InttMap::Offline_s, InttMap::OfflineWildcardComparator> Set_t;
  Set_t* m_bad_channel_set{nullptr};

  ClassDefOverride(InttBadChannelMapv1, 1)
};

#endif  // INTT_BAD_CHANNEL_MAPv1_H
