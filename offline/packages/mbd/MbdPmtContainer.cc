#include "MbdPmtContainer.h"
#include "MbdReturnCodes.h"
#include "MbdPmtHit.h"

#include <phool/phool.h>

#include <iostream>
#include <iomanip>
#include <cmath>

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

void MbdPmtContainer::Print(Option_t * /*option*/) const
{
  Short_t npmt = get_npmt();
  std::cout << "MBDPMT values:" << std::endl;
  for (Short_t ipmt=0; ipmt<npmt; ipmt++)
  {
    Float_t q = get_pmt(ipmt)->get_q();
    Float_t tt = get_pmt(ipmt)->get_tt();
    Float_t tq = get_pmt(ipmt)->get_tq();

    if ( std::fabs(tt)<100. )
    {
      std::cout << std::setw(8) << ipmt << std::setw(12) << q << "\t" << std::setw(12) << tt <<"\t" << std::setw(12) << tq << std::endl;
    }
  }
}

void MbdPmtContainer::virtual_warning(const std::string& funcsname) 
{
  std::cout << "MbdPmtContainer::" << funcsname << " is virtual, doing nothing" << std::endl;
  return;
}
