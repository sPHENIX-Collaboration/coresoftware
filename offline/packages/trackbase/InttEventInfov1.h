#ifndef INTT_EVENT_INFO_v1_H
#define INTT_EVENT_INFO_v1_H

#include "InttEventInfo.h"
#include "TrkrDefs.h"

#include <cstdint>
#include <cstddef>

#include <iostream>
#include <string>
#include <map>

#include <phool/PHObject.h>

class InttEventInfov1 : public InttEventInfo 
{
 public:
  InttEventInfov1();
  ~InttEventInfov1() override;

  void identify(std::ostream &os = std::cout) const override;
  void Reset() override;

  uint64_t get_bco_full() const override;
  void set_bco_full(uint64_t const&) override;

 protected:
  uint64_t bco_full;

 private:
  ClassDefOverride(InttEventInfov1, 1)
};

#endif//INTT_EVENT_INFO_v1_H
