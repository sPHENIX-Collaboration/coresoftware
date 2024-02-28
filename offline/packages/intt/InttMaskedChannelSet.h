#ifndef INTT_MASKED_CHANNEL_SET_H
#define INTT_MASKED_CHANNEL_SET_H

#include "InttMap.h"

#include <phool/PHObject.h>

#include <iostream>
#include <string>

class CDBTTree;

class InttMaskedChannelSet : public PHObject
{
 public:
  InttMaskedChannelSet() = default;
  ~InttMaskedChannelSet() override = default;

  virtual void identify(std::ostream& = std::cout) const override;
  virtual std::size_t size() const;

  int LoadFromFile(std::string const& = "InttMaskedChannelSet.root");
  int LoadFromCDB(std::string const& = "InttMaskedChannelSet");

  bool IsDeadChannel(int const&, int const&, int const&, int const&, int const&) const;
  virtual bool IsDeadChannel(InttMap::Offline_s const&) const;

 protected:
  virtual int v_LoadFromCDBTTree(CDBTTree&);

 private:
  ClassDefOverride(InttMaskedChannelSet, 1)
};

#endif  // INTT_MASKED_CHANNEL_SET_H
