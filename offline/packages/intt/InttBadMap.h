#ifndef INTT_BAD_MAP_H
#define INTT_BAD_MAP_H

#include "InttLoadable.h"
#include "InttMap.h"

#include <iostream>
#include <set>
#include <string>

class CDBTTree;

class InttBadMap : public InttLoadable
{
 public:
  InttBadMap() = default;
  virtual ~InttBadMap() override;

  bool IsBad(InttMap::Online_s const&) const;
  bool IsBad(InttMap::Offline_s const&) const;
  bool IsBad(InttMap::RawData_s const&) const;

  using InttLoadable::LoadFromFile;
  using InttLoadable::LoadFromCDB;

  int LoadFromFile() override { return LoadFromFile("InttBadMap.root"); }
  int LoadFromCDB() override { return LoadFromCDB("InttBadMap"); }

 protected:
  int LoadFromCDBTTree(CDBTTree&) override;

 private:
  typedef std::set<InttMap::RawData_s, InttMap::RawDataWildcardComparator> Set_t;
  Set_t* m_bad_channel_set{nullptr};
};

#endif  // INTT_BAD_MAP_H
