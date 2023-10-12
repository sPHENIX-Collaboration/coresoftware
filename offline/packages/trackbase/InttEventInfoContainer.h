#ifndef INTT_EVENT_INFO_CONTAINER_H
#define INTT_EVENT_INFO_CONTAINER_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <string>
#include <map>

class InttEventInfoContainer : public PHObject
{
public:
  struct InttEventInfo_s
  {
    Long64_t bco;
  };

  using KEY_t = Long64_t;
  using VAL_t = InttEventInfo_s;
  using Map = std::map<KEY_t, VAL_t>;
  using Iterator = Map::iterator;
  using ConstIterator = Map::const_iterator;

  InttEventInfoContainer();
  virtual ~InttEventInfoContainer() override;

  virtual void identify(std::ostream &os = std::cout) const override;
  virtual void Reset() override;

  virtual void AddInfo(KEY_t const&, VAL_t const&);
  virtual VAL_t& GetInfo(KEY_t const&);

 private:
  ClassDefOverride(InttEventInfoContainer, 1)
};

#endif//INTT_EVENT_INFO_CONTAINER_H
