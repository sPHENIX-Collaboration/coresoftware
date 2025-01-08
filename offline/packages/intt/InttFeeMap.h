#ifndef INTT_FEE_MAP_H
#define INTT_FEE_MAP_H

#include "InttMap.h"

#include <phool/PHObject.h>

#include <cstddef>  // for size_t
#include <iostream>
#include <string>

class CDBTTree;

class InttFeeMap : public PHObject
{
 public:
  InttFeeMap() = default;
  ~InttFeeMap() override = default;

  int LoadFromFile(std::string const& = "InttFeeMap.root");
  int LoadFromCDB(std::string const& = "InttFeeMap");

  virtual void identify(std::ostream& = std::cout) const override;

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
  ClassDefOverride(InttFeeMap, 1)
};

#endif  // INTT_FEE_MAP_H
