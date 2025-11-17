#include "MbdRawContainer.h"
#include "MbdReturnCodes.h"
#include "MbdRawHit.h"

#include <phool/phool.h>

#include <iostream>
#include <iomanip>
#include <cmath>

void MbdRawContainer::identify(std::ostream& os) const
{
  os << "virtual MbdRawContainer object" << std::endl;
  return;
}

void MbdRawContainer::Reset()
{
  std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
  return;
}

int MbdRawContainer::isValid() const
{
  virtual_warning("isValid()");
  return 0;
}

void MbdRawContainer::set_clocks(const Int_t /*evt*/, const UShort_t /*iclk*/, const UShort_t /*ifemclk*/)
{
  virtual_warning("set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk");
  return;
}

Int_t MbdRawContainer::get_evt() const
{
  virtual_warning("get_evt()");
  return 0;
}

UShort_t MbdRawContainer::get_clock() const
{
  virtual_warning("get_clock()");
  return 0;
}

UShort_t MbdRawContainer::get_femclock() const
{
  virtual_warning("get_femclock(const int iarm)");
  return 0;
}

void MbdRawContainer::set_npmt(const Short_t /*ival*/)
{
  virtual_warning("set_npmt(const Short_t ival)");
  return;
}

Short_t MbdRawContainer::get_npmt() const
{
  virtual_warning("get_npmt()");
  return MbdReturnCodes::MBD_INVALID_SHORT;
}

MbdRawHit *MbdRawContainer::get_pmt(const int /*iPmt*/) const
{
  virtual_warning("get_pmt(const short iPmt)");
  return nullptr;
}

void MbdRawContainer::Print(Option_t * /*option*/) const
{
  Short_t npmt = get_npmt();
  std::cout << "MbdRaw, evt " << get_evt() << std::endl;
  std::cout << "clk\t" << std::setw(12) << get_clock() << std::setw(12) << get_femclock() << std::endl;
  for (Short_t ipmt=0; ipmt<npmt; ipmt++)
  {
    Float_t adc = get_pmt(ipmt)->get_adc();
    Float_t ttdc = get_pmt(ipmt)->get_ttdc();
    Float_t qtdc = get_pmt(ipmt)->get_qtdc();

    std::cout << std::setw(8) << ipmt << std::setw(12) << adc << "\t" << std::setw(12) << ttdc <<"\t" << std::setw(12) << qtdc << std::endl;
  }
}

void MbdRawContainer::virtual_warning(const std::string& funcsname) 
{
  std::cout << "MbdRawContainer::" << funcsname << " is virtual, doing nothing" << std::endl;
  return;
}
