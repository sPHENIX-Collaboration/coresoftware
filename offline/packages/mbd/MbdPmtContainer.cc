#include "MbdPmtContainer.h"
#include "MbdReturnCodes.h"

#include <phool/phool.h>

#include <iostream>

void MbdPmtContainer::identify(std::ostream& os) const
{
  os << "virtual MbdPmtContainer object" << std::endl;
  return;
}

void MbdPmtContainer::Reset()
{
  std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
  return;
}

int MbdPmtContainer::isValid() const
{
  virtual_warning("isValid()");
  return 0;
}

void MbdPmtContainer::set_npmt(const Short_t /*ival*/)
{
  virtual_warning("set_npmt(const Short_t ival)");
  return;
}

Short_t MbdPmtContainer::get_npmt() const
{
  virtual_warning("get_npmt()");
  return MbdReturnCodes::MBD_INVALID_SHORT;
}

MbdPmtHit *MbdPmtContainer::get_pmt(const int /*iPmt*/) const
{
  virtual_warning("get_pmt(const short iPmt)");
  return nullptr;
}

void MbdPmtContainer::virtual_warning(const std::string& funcsname) const
{
  std::cout << "MbdPmtContainer::" << funcsname << " is virtual, doing nothing" << std::endl;
  return;
}
