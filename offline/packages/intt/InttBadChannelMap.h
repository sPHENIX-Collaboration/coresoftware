#ifndef INTT_BAD_CHANNEL_MAP_H
#define INTT_BAD_CHANNEL_MAP_H

#include "InttMap.h"

#include <iostream>
#include <string>
#include <set>

class CDBTTree;

class InttBadChannelMap
{
 public:
  InttBadChannelMap() = default;
  virtual ~InttBadChannelMap();

  int LoadFromFile(std::string const& = "InttBadChannelMap.root");
  int LoadFromCDB(std::string const& = "InttBadChannelMap");

  virtual void identify(std::ostream& = std::cout) const;
  virtual std::size_t size() const;

  virtual bool IsBad(InttMap::Online_s const&) const;
  virtual bool IsBad(InttMap::Offline_s const&) const;
  virtual bool IsBad(InttMap::RawData_s const&) const;

 protected:
  virtual int v_LoadFromCDBTTree(CDBTTree&);

 private:
  typedef std::set<InttMap::Offline_s, InttMap::OfflineWildcardComparator> Set_t;
  Set_t* m_bad_channel_set{nullptr};

};

#endif  // INTT_BAD_CHANNEL_MAP_H
