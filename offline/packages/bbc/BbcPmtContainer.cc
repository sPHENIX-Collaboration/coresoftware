#include "BbcPmtContainer.h"
#include "BbcReturnCodes.h"

#include <phool/phool.h>

#include <iostream>

void BbcPmtContainer::identify(std::ostream& os) const
{
  os << "virtual BbcPmtContainer object" << std::endl;
  return;
}

void BbcPmtContainer::Reset()
{
  std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
  return;
}

int BbcPmtContainer::isValid() const
{
  virtual_warning("isValid()");
  return 0;
}

void BbcPmtContainer::set_npmt(const short /*ival*/)
{
  virtual_warning("set_npmt(const short ival)");
  return;
}

short BbcPmtContainer::get_npmt() const
{
  virtual_warning("get_npmt()");
  return BbcReturnCodes::BBC_INVALID_SHORT;
}

short BbcPmtContainer::get_pmt(const int /*iPmt*/) const
{
  virtual_warning("get_pmt(const short iPmt)");
  return BbcReturnCodes::BBC_INVALID_SHORT;
}

float BbcPmtContainer::get_adc(const int /*iPmt*/) const
{
  virtual_warning("get_Adc(const short iPmt)");
  return BbcReturnCodes::BBC_INVALID_FLOAT;
}

float BbcPmtContainer::get_tdc0(const int /*iPmt*/) const
{
  virtual_warning("get_Tdc0(const short iPmt)");
  return BbcReturnCodes::BBC_INVALID_FLOAT;
}

float BbcPmtContainer::get_tdc1(const int /*iPmt*/) const
{
  virtual_warning("get_Tdc1(const short iPmt)");
  return BbcReturnCodes::BBC_INVALID_FLOAT;
}

void BbcPmtContainer::AddBbcPmt(const short /*ipmt*/, const float /*adc*/, const float /*tdc0*/, const float /*tdc1*/)
{
  virtual_warning("AddBbcPmtHit(const short ipmt, const short adc, const float tdc0, const float tdc1)");
  return;
}

void BbcPmtContainer::virtual_warning(const std::string& funcsname) const
{
  std::cout << "BbcPmtContainer::" << funcsname << " is virtual, doing nothing" << std::endl;
  return;
}
