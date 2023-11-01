#ifndef INTT_EVENT_INFO_CONTAINER_H
#define INTT_EVENT_INFO_CONTAINER_H

#include "TrkrDefs.h"

#include <cstdint>
#include <cstddef>

#include <iostream>
#include <string>
#include <map>

#include <phool/PHObject.h>

class InttEventInfoContainer : public PHObject
{
 public:
  InttEventInfoContainer();
  virtual ~InttEventInfoContainer();

  virtual void identify(std::ostream &os = std::cout) const;
  virtual void Reset();

  virtual uint64_t get_bco_full() const;
  virtual void set_bco_full(uint64_t const&);

 protected:

 private:

};

#endif//INTT_EVENT_INFO_CONTAINER_H
