#ifndef INTT_FEE_MAP_H
#define INTT_FEE_MAP_H

#include "InttMap.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <string>
#include <map>

class CDBTTree;

class InttFeeMap
{
 public:
  InttFeeMap() = default;
  virtual ~InttFeeMap();

  int LoadFromFile(std::string const& = "InttFeeMap.root");
  int LoadFromCDB(std::string const& = "InttFeeMap");

  virtual void identify(std::ostream& = std::cout) const;

  virtual int Convert(InttMap::Online_s&, InttMap::Offline_s const&) const;
  virtual int Convert(InttMap::Offline_s&, InttMap::Online_s const&) const;

  virtual int Convert(InttMap::RawData_s&, InttMap::Online_s const&) const;
  virtual int Convert(InttMap::Online_s&, InttMap::RawData_s const&) const;

  virtual int Convert(InttMap::RawData_s&, InttMap::Offline_s const&) const;
  virtual int Convert(InttMap::Offline_s&, InttMap::RawData_s const&) const;

  InttMap::Online_s ToOnline(InttMap::Offline_s const&) const;
  InttMap::Online_s ToOnline(InttMap::RawData_s const&) const;

  InttMap::RawData_s ToRawData(InttMap::Online_s const&) const;
  InttMap::RawData_s ToRawData(InttMap::Offline_s const&) const;

  InttMap::Offline_s ToOffline(InttMap::Online_s const&) const;
  InttMap::Offline_s ToOffline(InttMap::RawData_s const&) const;

 protected:
  virtual int v_LoadFromCDBTTree(CDBTTree&);

 private:
  std::map<InttMap::Online_s, InttMap::RawData_s, InttMap::OnlineWildcardComparator>* m_onl_to_raw{nullptr};
  std::map<InttMap::RawData_s, InttMap::Online_s, InttMap::RawDataWildcardComparator>* m_raw_to_onl{nullptr};
};

#endif  // INTT_FEE_MAP_H
