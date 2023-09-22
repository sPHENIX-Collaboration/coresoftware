#ifndef TRACKBASE_EVT_HEADER_V1_H
#define TRACKBASE_EVT_HEADER_V1_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <cstdint>

#include <iostream>
#include <map>

class EvtHeaderv1 : public EvtHeader
{
 public:
  EvtHeaderv1() = default;
  ~EvtHeaderv1() override {}
  void identify(std::ostream& = std::cout) const override;
  void Reset() override;

  void AddBCO(TrkrDefs::TrkrId const&, FEE_t const&, BCO_t const&) override;
  ConstIterator GetBCO(TrkrDefs::TrkrId const&, FEE_t const&) const override;

 protected:
  Map bco_map;

  ClassDefOverride(EvtHeaderv1, 1);
};

#endif//TRACKBASE_EVT_HEADER_V1_H
