#ifndef INTT_BAD_CHANNEL_MAP_H
#define INTT_BAD_CHANNEL_MAP_H

#include "InttMap.h"

#include <phool/PHObject.h>
#include <iostream>
#include <string>

class CDBTTree;

class InttBadChannelMap : public PHObject
{
 public:
  InttBadChannelMap() = default;
  ~InttBadChannelMap() override = default;

  int LoadFromFile(std::string const& = "InttBadChannelMap.root");
  int LoadFromCDB(std::string const& = "InttBadChannelMap");

  virtual void identify(std::ostream& = std::cout) const override;
  virtual std::size_t size() const;

  virtual bool IsBad(InttMap::Online_s const&) const;
  virtual bool IsBad(InttMap::Offline_s const&) const;
  virtual bool IsBad(InttMap::RawData_s const&) const;

 protected:
  virtual int v_LoadFromCDBTTree(CDBTTree&);

 private:
  ClassDefOverride(InttBadChannelMap, 1)
};

#endif  // INTT_BAD_CHANNEL_MAP_H
