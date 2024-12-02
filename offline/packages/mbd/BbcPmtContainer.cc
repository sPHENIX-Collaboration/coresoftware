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

void BbcPmtContainer::set_npmt(const Short_t /*ival*/)
{
  virtual_warning("set_npmt(const Short_t ival)");
  return;
}

Short_t BbcPmtContainer::get_npmt() const
{
  virtual_warning("get_npmt()");
  return BbcReturnCodes::BBC_INVALID_SHORT;
}

BbcPmtHit *BbcPmtContainer::get_pmt(const int /*iPmt*/) const
{
  virtual_warning("get_pmt(const short iPmt)");
  return nullptr;
}

void BbcPmtContainer::virtual_warning(const std::string& funcsname) const
{
  std::cout << "BbcPmtContainer::" << funcsname << " is virtual, doing nothing" << std::endl;
  return;
}
