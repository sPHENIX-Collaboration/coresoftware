#ifndef __BBCPMTHIT_H__
#define __BBCPMTHIT_H__

#include <iostream>
#include <phool/phool.h>
#include <TObject.h>

class BbcPmtHit : public TObject
{
public:

  BbcPmtHit() {}
  BbcPmtHit(const Short_t /*pmt*/, const Float_t /*adc*/, const Float_t /*tdc0*/, const Float_t /*tdc1*/) {}
  virtual ~BbcPmtHit() { }

  virtual Short_t get_pmt()  const { PHOOL_VIRTUAL_WARNING; return -9999; }
  virtual Float_t get_adc()  const { PHOOL_VIRTUAL_WARNING; return -9999; }
  virtual Float_t get_tdc0() const { PHOOL_VIRTUAL_WARNING; return -9999; }
  virtual Float_t get_tdc1() const { PHOOL_VIRTUAL_WARNING; return -9999; }

  void identify(std::ostream& os = std::cout) const;

private:
  ClassDefOverride(BbcPmtHit,1)
};

#endif
