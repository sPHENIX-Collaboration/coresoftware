#ifndef INTT_BAD_CHANNEL_MAP_H
#define INTT_BAD_CHANNEL_MAP_H

#include "InttLoadable.h"
#include "InttMap.h"

#include <set>
#include <string>

class CDBTTree;

class InttBadChannelMap : public InttLoadable
{
 public:
  InttBadChannelMap() = default;
  ~InttBadChannelMap() override = default;

  virtual bool IsBad(InttMap::Offline_s const&) const;

 protected:
  int LoadFromCdbTTree(CDBTTree&) override;

 private:
  typedef std::set<InttMap::Offline_s, InttMap::OfflineWildcardComparator> Set_t;
  Set_t m_bad_channel_set;
};

#endif  // INTT_BAD_CHANNEL_MAP_H
