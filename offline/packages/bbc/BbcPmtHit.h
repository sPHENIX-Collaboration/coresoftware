// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCPMTHIT_H
#define BBC_BBCPMTHIT_H

#include "BbcReturnCodes.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <iostream>

class BbcPmtHit : public PHObject
{
 public:
  BbcPmtHit() {}
  virtual ~BbcPmtHit() override = default;

  virtual Short_t get_pmt() const
  {
    PHOOL_VIRTUAL_WARNING;
    return -9999;
  }

  virtual Float_t get_q() const
  {
    PHOOL_VIRTUAL_WARNING;
    return BbcReturnCodes::BBC_INVALID_FLOAT;
  }

  virtual Float_t get_time() const
  {
    PHOOL_VIRTUAL_WARNING;
    return BbcReturnCodes::BBC_INVALID_FLOAT;
  }

  virtual Float_t get_tt() const
  {
    PHOOL_VIRTUAL_WARNING;
    return BbcReturnCodes::BBC_INVALID_FLOAT;
  }

  virtual Float_t get_tq() const
  {
    PHOOL_VIRTUAL_WARNING;
    return BbcReturnCodes::BBC_INVALID_FLOAT;
  }

  virtual void set_pmt(const Short_t /*pmt*/, const Float_t /*q*/, const Float_t /*tt*/, const Float_t /*tq*/)
  {
    PHOOL_VIRTUAL_WARNING;
  }

  virtual void identify(std::ostream& os = std::cout) const override;

  virtual int isValid() const override { return 0; }

 private:
  ClassDefOverride(BbcPmtHit, 1)
};

#endif
