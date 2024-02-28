#ifndef INTT_MASKED_CHANNEL_SETv1_H
#define INTT_MASKED_CHANNEL_SETv1_H

#include "InttMap.h"
#include "InttMaskedChannelSet.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <set>

class CDBTTree;

class InttMaskedChannelSetv1 : public InttMaskedChannelSet
{
 public:
  InttMaskedChannelSetv1() = default;
  ~InttMaskedChannelSetv1() override;

  void identify(std::ostream& = std::cout) const override;
  std::size_t size() const override;

  bool IsDeadChannel(InttMap::Offline_s const&) const override;

 protected:
  int v_LoadFromCDBTTree(CDBTTree&) override;

 private:
  typedef std::set<InttMap::Offline_s, InttMap::OfflineWildcardComparator> Set_t;
  Set_t* m_HotChannelSet{nullptr};

  ClassDefOverride(InttMaskedChannelSetv1, 1)
};

#endif  // INTT_MASKED_CHANNEL_SETv1_H
