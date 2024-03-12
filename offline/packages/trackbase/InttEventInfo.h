#ifndef INTT_EVENT_INFO_H
#define INTT_EVENT_INFO_H

#include "TrkrDefs.h"

#include <cstddef>
#include <cstdint>

#include <iostream>
#include <map>
#include <string>

#include <phool/PHObject.h>

class InttEventInfo : public PHObject
{
 public:
  InttEventInfo() = default;
  virtual ~InttEventInfo() = default;

  virtual void identify(std::ostream &os = std::cout) const override;
  virtual void Reset() override;

  virtual uint64_t get_bco_full() const;
  virtual void set_bco_full(uint64_t const &);

 protected:
  ClassDefOverride(InttEventInfo, 1);
};

#endif  // INTT_EVENT_INFO_H
