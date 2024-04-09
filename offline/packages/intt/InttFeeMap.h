#ifndef INTT_FEE_MAP_H
#define INTT_FEE_MAP_H

#include "InttLoadable.h"
#include "InttMap.h"

#include <map>
#include <string>

class CDBTTree;

class InttFeeMap : public InttLoadable
{
 public:
  InttFeeMap() = default;
  ~InttFeeMap() override;

  int Convert(InttMap::Online_s&, InttMap::Offline_s const&) const;
  int Convert(InttMap::Offline_s&, InttMap::Online_s const&) const;

  int Convert(InttMap::RawData_s&, InttMap::Online_s const&) const;
  int Convert(InttMap::Online_s&, InttMap::RawData_s const&) const;

  int Convert(InttMap::RawData_s&, InttMap::Offline_s const&) const;
  int Convert(InttMap::Offline_s&, InttMap::RawData_s const&) const;

  InttMap::Online_s ToOnline(InttMap::Offline_s const&) const;
  InttMap::Online_s ToOnline(InttMap::RawData_s const&) const;

  InttMap::RawData_s ToRawData(InttMap::Online_s const&) const;
  InttMap::RawData_s ToRawData(InttMap::Offline_s const&) const;

  InttMap::Offline_s ToOffline(InttMap::Online_s const&) const;
  InttMap::Offline_s ToOffline(InttMap::RawData_s const&) const;

  using InttLoadable::LoadFromFile;
  using InttLoadable::LoadFromCDB;

  int LoadFromFile() override { return LoadFromFile("InttBadChannelMap.root"); }
  int LoadFromCDB() override { return LoadFromCDB("InttBadChannelMap"); }

 protected:
  int LoadFromCDBTTree(CDBTTree&) override;

 private:
  typedef std::map<InttMap::Online_s, InttMap::RawData_s, InttMap::OnlineWildcardComparator> onl_to_raw_map_t;
  typedef std::map<InttMap::RawData_s, InttMap::Online_s, InttMap::RawDataWildcardComparator> raw_to_onl_map_t;

  onl_to_raw_map_t* m_onl_to_raw{nullptr};
  raw_to_onl_map_t* m_raw_to_onl{nullptr};
};

#endif  // INTT_FEE_MAP_H
